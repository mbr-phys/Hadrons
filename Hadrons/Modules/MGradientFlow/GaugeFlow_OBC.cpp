/*
 * GaugeFlow_OBC.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Matthew Black    <matthewkblack@protonmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 */

/*  END LEGAL */
#include <Hadrons/Modules/MGradientFlow/GaugeFlow_OBC.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGradientFlow;

template class Grid::Hadrons::MGradientFlow::TGaugeFlow_OBC<GIMPL, OBC::WilsonAction<GIMPL>>;
template class Grid::Hadrons::MGradientFlow::TGaugeFlow_OBC<GIMPL, SymanzikGaugeAction<GIMPL>>;
template class Grid::Hadrons::MGradientFlow::TGaugeFlow_OBC<GIMPL, PlaqPlusRectangleAction<GIMPL>>;
