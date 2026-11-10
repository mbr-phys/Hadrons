/*
 * MomentumProjectMany.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Matthew Black <matthewkblack@protonmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MomentumProjectMany_hpp_
#define Hadrons_MomentumProjectMany_hpp_

#include <Hadrons/Global.hpp>
#include <Grid/algorithms/blas/MomentumProject.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         MomentumProjectMany                                *
 ******************************************************************************
 * Grid's MomentumProject, extended to project several fields at once.
 *
 * ProjectMany packs the fields as extra components of a single BLAS matrix,
 *   V[xyz][t][field][word]  (Fortran column-major V_(word,field,t),xyz),
 * so that one GEMM against the phase matrix M_xyz,mom of the base class, and
 * one global sum, project all of them onto every momentum. Project() handles
 * one field per call, with a GEMM of only words*nt rows (nt local); packing
 * nField fields multiplies that dimension by nField and replaces nField
 * imports, GEMMs and global sums by one of each.
 *
 * The projected fields may have any type with the scalar type of the base
 * class (e.g. LatticeComplex traces, one per gamma structure, or colour
 * matrices). The phases are those imported by the base class, so the
 * momentum list and phase sign are those of the projector module.
 *
 * Output: out[field][t + nt*mom], t global, as in MomentumProject::Project.
 */
template <class Field, class ComplexField>
class MomentumProjectMany: public MomentumProject<Field, ComplexField>
{
public:
    typedef MomentumProject<Field, ComplexField> Base;
    typedef typename Base::scalar                scalar;
public:
    template <class ProjField>
    void ProjectMany(const std::vector<ProjField> &fields,
                     std::vector<std::vector<typename ProjField::scalar_object>> &out);
private:
    template <class ProjField>
    void importVectors(const std::vector<ProjField> &fields);
private:
    deviceVector<scalar> manyV_, manyP_;
};

/******************************************************************************
 *                 MomentumProjectMany implementation                         *
 ******************************************************************************/
// V[xyz][t][field][word], the layout of MomentumProject::ImportVector with the
// fields as an extra index between t and word
template <class Field, class ComplexField>
template <class ProjField>
void MomentumProjectMany<Field, ComplexField>::importVectors(const std::vector<ProjField> &fields)
{
    typedef typename ProjField::vector_object vobj;
    typedef typename ProjField::scalar_object sobj;

    GridBase         *grid        = this->grid;
    const int        nd           = grid->_ndimension;
    const Coordinate rdimensions  = grid->_rdimensions;
    const Coordinate ldims        = grid->LocalDimensions();
    const Coordinate simd         = grid->_simd_layout;
    const int64_t    osites       = grid->oSites();
    const int        Nsimd        = vobj::Nsimd();
    const int64_t    nwords       = sizeof(sobj)/sizeof(scalar);
    const int64_t    nField       = fields.size();
    const int64_t    Nt           = this->nt;

    manyV_.resize(this->nxyz*Nt*nField*nwords);
    auto blasData_p = &manyV_[0];
    for (int64_t f = 0; f < nField; ++f)
    {
        GRID_ASSERT(fields[f].Grid() == grid);
        autoView(Data, fields[f], AcceleratorRead);
        accelerator_for(sf, osites, Nsimd, {
#ifdef GRID_SIMT
        {
            int lane = acceleratorSIMTlane(Nsimd); // buffer lane
#else
        for (int lane = 0; lane < Nsimd; lane++) {
#endif
            Coordinate lcoor(nd, 0);
            Coordinate icoor(nd);
            Coordinate ocoor(nd);

            Lexicographic::CoorFromIndex(icoor, lane, simd);
            Lexicographic::CoorFromIndex(ocoor, sf, rdimensions);
            for (int d = 0; d < nd; d++)
            {
                lcoor[d] = rdimensions[d]*icoor[d] + ocoor[d];
            }
            uint64_t l_t = lcoor[nd - 1];

            int64_t    l_xyz    = 0;
            Coordinate xyz_coor = lcoor;
            xyz_coor[nd - 1] = 0;
            Lexicographic::IndexFromCoor(xyz_coor, l_xyz, ldims);

            sobj   data       = extractLane(lane, Data[sf]);
            scalar *data_words = (scalar *) &data;
            for (int64_t w = 0; w < nwords; w++)
            {
                uint64_t idx = w + nwords*(f + nField*(l_t + Nt*l_xyz));
                blasData_p[idx] = data_words[w];
            }
        }
        });
    }
}

template <class Field, class ComplexField>
template <class ProjField>
void MomentumProjectMany<Field, ComplexField>::ProjectMany(
    const std::vector<ProjField> &fields,
    std::vector<std::vector<typename ProjField::scalar_object>> &out)
{
    typedef typename ProjField::scalar_object sobj;
    static_assert(std::is_same<typename ProjField::scalar_type, scalar>::value,
                  "projected fields must have the scalar type of the projector");

    GridBase      *grid  = this->grid;
    const int     nd     = grid->Nd();
    const int64_t nwords = sizeof(sobj)/sizeof(scalar);
    const int64_t nField = fields.size();
    const int64_t lt     = this->nt;
    const int64_t gt     = grid->GlobalDimensions()[nd - 1];
    const int64_t st     = grid->LocalStarts()[nd - 1];
    const int64_t nmom   = this->nmom;
    const int64_t rows   = nwords*nField*lt;
    sobj          zero;

    zeroit(zero);
    out.assign(nField, std::vector<sobj>(gt*nmom, zero));
    if (nField == 0)
    {
        return;
    }
    importVectors(fields);
    manyP_.resize(rows*nmom);

    deviceVector<scalar *> Vd(1);
    deviceVector<scalar *> Md(1);
    deviceVector<scalar *> Pd(1);
    scalar *Vh = &manyV_[0];
    scalar *Mh = &this->BLAS_M[0];
    scalar *Ph = &manyP_[0];
    acceleratorPut(Vd[0], Vh);
    acceleratorPut(Md[0], Mh);
    acceleratorPut(Pd[0], Ph);

    // P_(word,field,t),mom = V_(word,field,t),xyz M_xyz,mom
    GridBLAS BLAS;
    BLAS.gemmBatched(GridBLAS_OP_N, GridBLAS_OP_N,
                     rows, nmom, this->nxyz,
                     scalar(1.0),
                     Vd,
                     Md,
                     scalar(0.0),
                     Pd);
    BLAS.synchronise();

    // local P[mom][t][field][word] -> global [field][mom][t][word], zero
    // outside this rank's time slices, then one sum over the ranks
    std::vector<scalar> planes(rows*nmom);
    std::vector<scalar> gdata(nField*nmom*gt*nwords, scalar(0.0));

    acceleratorCopyFromDevice(&manyP_[0], &planes[0], planes.size()*sizeof(scalar));
    for (int64_t m = 0; m < nmom; ++m)
    for (int64_t t = 0; t < lt; ++t)
    for (int64_t f = 0; f < nField; ++f)
    for (int64_t w = 0; w < nwords; ++w)
    {
        gdata[w + nwords*((t + st) + gt*(m + nmom*f))] =
            planes[w + nwords*(f + nField*(t + lt*m))];
    }
    grid->GlobalSumVector(&gdata[0], gdata.size());
    for (int64_t f = 0; f < nField; ++f)
    for (int64_t i = 0; i < nmom*gt; ++i)
    {
        scalar *o = (scalar *) &out[f][i];

        for (int64_t w = 0; w < nwords; ++w)
        {
            o[w] = gdata[w + nwords*(i + nmom*gt*f)];
        }
    }
}

END_HADRONS_NAMESPACE

#endif // Hadrons_MomentumProjectMany_hpp_
