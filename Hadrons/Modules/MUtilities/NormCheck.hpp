/*
 * NormCheck.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2026
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 */
#ifndef Hadrons_MUtilities_NormCheck_hpp_
#define Hadrons_MUtilities_NormCheck_hpp_

#include <algorithm>
#include <cmath>

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MUtilities)

class NormCheckPar: Serializable
{
public:
    // reference is optional.  When supplied, the module also reports the
    // inner product and the relative norm of the difference.
    GRID_SERIALIZABLE_CLASS_MEMBERS(NormCheckPar,
                                    std::string, field,
                                    std::string, reference,
                                    std::string, output);
};

template <typename Field>
class TNormCheck: public Module<NormCheckPar>
{
public:
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        double,  norm2,
                                        double,  referenceNorm2,
                                        Complex, innerProduct,
                                        double,  differenceNorm2,
                                        double,  relativeDifference);
    };
public:
    TNormCheck(const std::string name);
    virtual ~TNormCheck(void) {};
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    virtual void setup(void);
    virtual void execute(void);
};

MODULE_REGISTER_TMP(NormCheckFermion, TNormCheck<FIMPL::FermionField>, MUtilities);
MODULE_REGISTER_TMP(NormCheckPropagator, TNormCheck<FIMPL::PropagatorField>, MUtilities);
MODULE_REGISTER_TMP(NormCheckComplex, TNormCheck<FIMPL::ComplexField>, MUtilities);
MODULE_REGISTER_TMP(NormCheckColourMatrix, TNormCheck<GIMPL::GaugeLinkField>, MUtilities);

template <typename Field>
TNormCheck<Field>::TNormCheck(const std::string name)
: Module<NormCheckPar>(name)
{}

template <typename Field>
std::vector<std::string> TNormCheck<Field>::getInput(void)
{
    std::vector<std::string> in = {par().field};

    if (!par().reference.empty())
    {
        in.push_back(par().reference);
    }
    return in;
}

template <typename Field>
std::vector<std::string> TNormCheck<Field>::getOutput(void)
{
    return {getName()};
}

template <typename Field>
std::vector<std::string> TNormCheck<Field>::getOutputFiles(void)
{
    if (par().output.empty())
    {
        return {};
    }
    return {resultFilename(par().output)};
}

template <typename Field>
void TNormCheck<Field>::setup(void)
{
    if (par().field.empty())
    {
        HADRONS_ERROR(Argument, "field must be specified");
    }
    envCreate(HadronsSerializable, getName(), 1, 0);
}

template <typename Field>
void TNormCheck<Field>::execute(void)
{
    const auto &field = envGet(Field, par().field);
    Result result{};

    result.norm2 = norm2(field);
    if (!par().reference.empty())
    {
        const auto &reference = envGet(Field, par().reference);
        result.referenceNorm2 = norm2(reference);
        result.innerProduct = innerProduct(field, reference);
        result.differenceNorm2 = result.norm2 + result.referenceNorm2
                               - 2.0*real(result.innerProduct);

        const double denominator = std::sqrt(result.norm2)
                                 + std::sqrt(result.referenceNorm2);
        result.relativeDifference = (denominator == 0.0)
                                  ? 0.0
                                  : std::sqrt(std::max(0.0, result.differenceNorm2))/denominator;
    }

    LOG(Message) << "NormCheck '" << par().field << "': norm2 = "
                 << result.norm2 << std::endl;
    if (!par().reference.empty())
    {
        LOG(Message) << "  reference '" << par().reference
                     << "': norm2 = " << result.referenceNorm2
                     << ", difference norm2 = " << result.differenceNorm2
                     << ", relative difference = " << result.relativeDifference
                     << std::endl;
    }

    saveResult(par().output, "norm_check", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_NormCheck_hpp_
