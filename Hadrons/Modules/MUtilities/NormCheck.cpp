/*
 * NormCheck.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2026
 */
#include <Hadrons/Modules/MUtilities/NormCheck.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MUtilities;

template class HADRONS_NAMESPACE::MUtilities::TNormCheck<FIMPL::FermionField>;
template class HADRONS_NAMESPACE::MUtilities::TNormCheck<FIMPL::PropagatorField>;
template class HADRONS_NAMESPACE::MUtilities::TNormCheck<FIMPL::ComplexField>;
template class HADRONS_NAMESPACE::MUtilities::TNormCheck<GIMPL::GaugeLinkField>;
