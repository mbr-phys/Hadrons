/*
 * Test_gradient_flow_boundary_comparison.cpp, part of Hadrons
 *
 * Copyright (C) 2026
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 */

#include <Hadrons/Application.hpp>
#include <Hadrons/GaugeBoundary.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());

    Application application;
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1000;
    globalPar.trajCounter.end   = 1000;
    globalPar.trajCounter.step  = 1;
    globalPar.runId             = "gradient-flow-bc-comparison";
    application.setPar(globalPar);

    application.createModule<MGauge::Random>("gauge");

    // (a) Existing periodic gauge-flow module.
    MGradientFlow::WilsonFlow::Par standardPar;
    standardPar.output         = "";
    standardPar.gauge          = "gauge";
    standardPar.steps          = 10;
    standardPar.step_size      = 0.01;
    standardPar.meas_interval  = 1;
    standardPar.save_history   = false;
    standardPar.save_rk_stages = false;
    application.createModule<MGradientFlow::WilsonFlow>("standard_periodic", standardPar);

    // (b) OBC-capable module operating with a periodic gauge boundary.
    MGradientFlow::WilsonFlow_OBC::Par obcPeriodicPar;
    obcPeriodicPar.output         = "";
    obcPeriodicPar.gauge          = "gauge";
    obcPeriodicPar.steps          = 10;
    obcPeriodicPar.step_size      = 0.01;
    obcPeriodicPar.meas_interval  = 1;
    obcPeriodicPar.save_history   = false;
    obcPeriodicPar.save_rk_stages = false;
    obcPeriodicPar.gaugeBC        = GaugeBoundary::periodic;
    obcPeriodicPar.cG             = 1.0;
    application.createModule<MGradientFlow::WilsonFlow_OBC>("obc_periodic", obcPeriodicPar);

    // (c) OBC-capable module operating with an open temporal gauge boundary.
    MGradientFlow::WilsonFlow_OBC::Par obcOpenPar;
    obcOpenPar.output         = "";
    obcOpenPar.gauge          = "gauge";
    obcOpenPar.steps          = 10;
    obcOpenPar.step_size      = 0.01;
    obcOpenPar.meas_interval  = 1;
    obcOpenPar.save_history   = false;
    obcOpenPar.save_rk_stages = false;
    obcOpenPar.gaugeBC        = GaugeBoundary::open;
    obcOpenPar.cG             = 1.0;
    application.createModule<MGradientFlow::WilsonFlow_OBC>("obc_open", obcOpenPar);

    MIO::WriteResultGroup::Par bundlePar;
    bundlePar.results = {"standard_periodic", "obc_periodic", "obc_open"};
    bundlePar.output  = "gradient_flow_boundary_comparison";
    application.createModule<MIO::WriteResultGroup>("write_comparison", bundlePar);

    application.saveParameterFile("GradientFlowBoundaryComparison.xml");
    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
