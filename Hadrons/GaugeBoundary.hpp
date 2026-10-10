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

BEGIN_HADRONS_NAMESPACE

namespace GaugeBoundary
{

inline constexpr const char *unknown  = "unknown";
inline constexpr const char *periodic = "periodic";
inline constexpr const char *open     = "open";

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

inline std::string gaugeTypeFromFermionBc(const int bc)
{
    return (bc == 0) ? open : periodic;
}

inline RealD openQcdPlaquetteWeight(const Coordinate &dimensions,
                                    const int timeDir)
{
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
                                   const int timeDir, const RealD cG)
{
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

class GaugeBoundaryInfo: public Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeBoundaryInfo,
                                    std::string, type,
                                    int, timeDir,
                                    double, cG,
                                    bool, userSpecified,
                                    bool, hasFlowInfo,
                                    std::string, flowAction,
                                    double, stepSize,
                                    double, flowTime,
                                    std::string, stage);

    GaugeBoundaryInfo(void)
    : type(GaugeBoundary::unknown)
    , timeDir(Nd - 1)
    , cG(1.0)
    , userSpecified(false)
    , hasFlowInfo(false)
    , stepSize(0.0)
    , flowTime(0.0)
    {}

    GaugeBoundaryInfo(const std::string &type_, const int timeDir_,
                      const double cG_, const bool userSpecified_,
                      const bool hasFlowInfo_ = false,
                      const std::string &flowAction_ = std::string(),
                      const double stepSize_ = 0.0,
                      const double flowTime_ = 0.0,
                      const std::string &stage_ = std::string())
    : type(type_)
    , timeDir(timeDir_)
    , cG(cG_)
    , userSpecified(userSpecified_)
    , hasFlowInfo(hasFlowInfo_)
    , flowAction(flowAction_)
    , stepSize(stepSize_)
    , flowTime(flowTime_)
    , stage(stage_)
    {}
};

namespace GaugeBoundary
{

inline bool isKnown(const GaugeBoundaryInfo &info)
{
    return info.type != unknown;
}

inline std::string compatibilityError(const GaugeBoundaryInfo &info,
                                      const std::string &expectedType,
                                      const int expectedTimeDir,
                                      const double expectedCG)
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
    if (expectedType == open) {
        if (info.timeDir != expectedTimeDir) {
            return "metadata declares open direction " +
                   std::to_string(info.timeDir) + " but module expects " +
                   std::to_string(expectedTimeDir);
        }
        const double scale = 1.0 + std::max(std::abs(info.cG),
                                            std::abs(expectedCG));
        if (!std::isfinite(info.cG) ||
            (std::abs(info.cG - expectedCG) > 1.0e-12 * scale)) {
            return "metadata declares cG=" + std::to_string(info.cG) +
                   " but module expects cG=" + std::to_string(expectedCG);
        }
    }

    return std::string();
}

} // namespace GaugeBoundary

END_HADRONS_NAMESPACE

#endif // Hadrons_GaugeBoundary_hpp_
