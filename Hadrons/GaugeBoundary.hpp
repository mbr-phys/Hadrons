/*
 * GaugeBoundary.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2026
 *
 * Author: Matthew Black <matthewkblack@protonmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
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
#ifndef Hadrons_GaugeBoundary_hpp_
#define Hadrons_GaugeBoundary_hpp_

#include <Hadrons/Serialization.hpp>
#include <algorithm>
#include <cmath>
#include <limits>

BEGIN_HADRONS_NAMESPACE

namespace GaugeBoundary
{

inline constexpr const char *unknown  = "unknown";
inline constexpr const char *periodic = "periodic";
inline constexpr const char *open     = "open";

// Open-boundary flow (openQCD wflow.c, Luescher & Schaefer CPC 184 (2013) 519,
// eqs. (B.2)-(B.4)): the force on the spatial links at x0 = 0 and x0 = T is
// multiplied by 1/w(x0) = 2. It is independent of cG; with flow cG = 1 the open
// flow equals the periodic flow of the lattice reflected across both
// boundaries (the orbifold construction mentioned there).
inline constexpr double openFlowBoundaryForceFactor = 2.0;

// cG value of an unknown ensemble (periodic gauge field or no declaration)
inline double unknownCG(void)
{
    return std::numeric_limits<double>::quiet_NaN();
}

inline bool closeCG(const double a, const double b)
{
    const double scale = 1.0 + std::max(std::abs(a), std::abs(b));
    return std::isfinite(a) && std::isfinite(b) &&
           (std::abs(a - b) <= 1.0e-12 * scale);
}

inline bool isValidType(const std::string &type)
{
    return (type == unknown) || (type == periodic) || (type == open);
}

inline bool isOpen(const std::string &type)
{
    return type == open;
}

inline std::string metadataName(const std::string &gaugeName)
{
    return gaugeName + "_boundary";
}

inline bool isMetadataName(const std::string &name)
{
    const std::string suffix = metadataName("");
    return (name.size() > suffix.size()) &&
           (name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0);
}

inline std::string gaugeTypeFromFermionBc(const int bc)
{
    return (bc == 0) ? open : periodic;
}

inline RealD openQcdPlaquetteWeight(const Coordinate &dimensions)
{
    constexpr int timeDir = Nd - 1;
    RealD transverseVolume = 1.0;
    for (int mu = 0; mu < static_cast<int>(dimensions.size()); ++mu) {
        if (mu != timeDir) {
            transverseVolume *= dimensions[mu];
        }
    }
    const RealD nd = dimensions.size();
    const RealD nt = dimensions[timeDir];
    const RealD nTemporalPlanes = nd - 1.0;
    const RealD nSpatialPlanes = (nd - 1.0) * (nd - 2.0) * 0.5;

    return transverseVolume * (nt - 1.0) *
           (nTemporalPlanes + nSpatialPlanes);
}

inline RealD actionPlaquetteWeight(const Coordinate &dimensions,
                                   const RealD cG)
{
    constexpr int timeDir = Nd - 1;
    RealD transverseVolume = 1.0;
    for (int mu = 0; mu < static_cast<int>(dimensions.size()); ++mu) {
        if (mu != timeDir) {
            transverseVolume *= dimensions[mu];
        }
    }
    const RealD nd = dimensions.size();
    const RealD nt = dimensions[timeDir];
    const RealD nTemporalPlanes = nd - 1.0;
    const RealD nSpatialPlanes = (nd - 1.0) * (nd - 2.0) * 0.5;

    return transverseVolume *
           (nTemporalPlanes * (nt - 1.0) +
            nSpatialPlanes * (nt - 2.0 + cG));
}

} // namespace GaugeBoundary

// cG is the boundary coefficient of the action the ensemble was generated with
// (NaN if not known) and is passed unchanged through flows. flowCG is the
// boundary coefficient of the flow action (meaningful when hasFlowInfo). Both
// use the openQCD weight cG/2 for spatial loops on the open boundary slices
// (CPC 184 (2013) 519, eq. (A.3)).
class GaugeBoundaryInfo: public Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeBoundaryInfo,
                                    std::string, type,
                                    double, cG,
                                    bool, userSpecified,
                                    bool, hasFlowInfo,
                                    std::string, flowAction,
                                    double, stepSize,
                                    double, flowTime,
                                    std::string, stage,
                                    double, flowCG);

    GaugeBoundaryInfo(void)
    : type(GaugeBoundary::unknown)
    , cG(1.0)
    , userSpecified(false)
    , hasFlowInfo(false)
    , stepSize(0.0)
    , flowTime(0.0)
    , flowCG(unknownCG())
    {}

    GaugeBoundaryInfo(const std::string &type_, const double cG_,
                      const bool userSpecified_,
                      const bool hasFlowInfo_ = false,
                      const std::string &flowAction_ = std::string(),
                      const double stepSize_ = 0.0,
                      const double flowTime_ = 0.0,
                      const std::string &stage_ = std::string(),
                      const double flowCG_ = unknownCG())
    : type(type_)
    , cG(cG_)
    , userSpecified(userSpecified_)
    , hasFlowInfo(hasFlowInfo_)
    , flowAction(flowAction_)
    , stepSize(stepSize_)
    , flowTime(flowTime_)
    , stage(stage_)
    , flowCG(flowCG_)
    {}
};

namespace GaugeBoundary
{

inline bool isKnown(const GaugeBoundaryInfo &info)
{
    return info.type != unknown;
}

// boundary type only: the ensemble cG and the flow cG are independent
inline std::string compatibilityError(const GaugeBoundaryInfo &info,
                                      const std::string &expectedType)
{
    if (!isValidType(info.type)) {
        return "metadata contains invalid gauge boundary type '" + info.type + "'";
    }
    if (!isValidType(expectedType) || (expectedType == unknown)) {
        return "internal expected gauge boundary type is invalid";
    }
    if (!isKnown(info)) {
        return std::string();
    }
    if (info.type != expectedType) {
        return "metadata declares gauge boundary '" + info.type +
               "' but module expects '" + expectedType + "'";
    }

    return std::string();
}

// flow provenance: a field produced by a flow must have been flowed with the
// expected action and, for open boundaries, the expected flow cG
inline std::string flowCompatibilityError(const GaugeBoundaryInfo &info,
                                          const std::string &expectedAction,
                                          const double expectedFlowCG)
{
    if (!info.hasFlowInfo) {
        return std::string();
    }
    if (info.flowAction != expectedAction) {
        return "metadata declares flow action '" + info.flowAction +
               "' but module expects '" + expectedAction + "'";
    }
    if ((info.type == open) && !closeCG(info.flowCG, expectedFlowCG)) {
        return "metadata declares flow cG=" + std::to_string(info.flowCG) +
               " but module expects flow cG=" + std::to_string(expectedFlowCG);
    }

    return std::string();
}

} // namespace GaugeBoundary

END_HADRONS_NAMESPACE

#endif // Hadrons_GaugeBoundary_hpp_
