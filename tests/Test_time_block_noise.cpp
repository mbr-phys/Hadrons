/*
 * Test_time_block_noise.cpp, part of Hadrons
 *
 * Copyright (C) 2026
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Checks TimeBlockNoise (MNoise::TimeBlockSpinColorDiagonal): every component
 * vanishes outside its block of time slices and not inside it, and the
 * components of a base noise sum to the undiluted spin-colour diagonal noise,
 * with no overlap (the norms add up). Use Nt >= 4.
 */
#include <Hadrons/Global.hpp>
#include <Hadrons/DilutedNoise.hpp>

using namespace Grid;
using namespace Hadrons;

typedef FIMPL                  FImpl;
typedef FImpl::PropagatorField PropagatorField;

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);

    GridCartesian   *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                GridDefaultSimd(Nd, vComplex::Nsimd()),
                                GridDefaultMpi());
    GridParallelRNG rng(grid);
    const int       nt    = grid->GlobalDimensions()[Tp];
    const int       nsrc  = 2;
    const double    tol   = 1.0e-12;
    bool            pass  = true;

    if (nt < 4)
    {
        std::cout << GridLogError << "Nt must be at least 4" << std::endl;
        Grid_finalize();

        return EXIT_FAILURE;
    }

    // three blocks of unequal length
    const std::vector<unsigned int> starts = {0, static_cast<unsigned int>(nt/4),
                                              static_cast<unsigned int>(nt/2 + 1)};
    TimeBlockNoise<FImpl> noise(grid, nsrc, starts);
    const unsigned int    nBlock = noise.nBlock();

    rng.SeedFixedIntegers({1, 2, 3, 4});
    noise.generateNoise(rng);
    if (noise.dilutionSize() != static_cast<int>(nsrc*nBlock))
    {
        std::cout << GridLogError << "FAIL: dilution size " << noise.dilutionSize()
                  << ", expected " << nsrc*nBlock << std::endl;
        pass = false;
    }

    PropagatorField       comp(grid), sum(grid), full(grid);
    LatticeComplex        n2(grid);
    std::vector<TComplex> slices;

    for (int r = 0; r < nsrc; ++r)
    {
        double normSum = 0.;

        sum = Zero();
        for (unsigned int b = 0; b < nBlock; ++b)
        {
            const int begin = noise.blockBegin(b), end = noise.blockEnd(b);

            comp     = noise.getProp(r*nBlock + b);
            sum      = sum + comp;
            normSum += norm2(comp);

            // support: non-zero exactly on [begin, end)
            n2 = localNorm2(comp);
            sliceSum(n2, slices, Tp);
            for (int t = 0; t < nt; ++t)
            {
                const double s      = real(TensorRemove(slices[t]));
                const bool   inside = (t >= begin) and (t < end);

                if (inside ? (s == 0.) : (s != 0.))
                {
                    std::cout << GridLogError << "FAIL: noise " << r << ", block " << b
                              << " [" << begin << ", " << end << "), slice " << t
                              << ": norm " << s << std::endl;
                    pass = false;
                }
            }
        }

        // the blocks partition the time extent
        full = 1.;
        full = full*noise.getNoise()[r];

        const double dSum  = norm2(sum - full)/norm2(full);
        const double dNorm = std::abs(normSum - norm2(full))/norm2(full);

        std::cout << GridLogMessage << "noise " << r << ": |sum of blocks - noise|^2/|noise|^2 = "
                  << dSum << ", relative norm defect = " << dNorm << std::endl;
        if (!(dSum < tol) or !(dNorm < tol))
        {
            std::cout << GridLogError << "FAIL: blocks of noise " << r
                      << " do not partition it" << std::endl;
            pass = false;
        }
    }

    std::cout << GridLogMessage << (pass ? "PASS" : "FAIL") << std::endl;
    Grid_finalize();

    return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
