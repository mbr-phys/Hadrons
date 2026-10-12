/*
 * RBPrecBiCGSTAB.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2026
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 */
#include <Hadrons/Modules/MSolver/RBPrecBiCGSTAB.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MSolver;

template class HADRONS_NAMESPACE::MSolver::TRBPrecBiCGSTAB<FIMPL>;
template class HADRONS_NAMESPACE::MSolver::TRBPrecBiCGSTAB<ZFIMPL>;
