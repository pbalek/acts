// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Propagator/MultiEigenStepperLoop.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/GaussianSumFitter.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"

#include <filesystem>

using namespace ActsExamples;

namespace {

using MultiStepper = Acts::MultiEigenStepperLoop<>;
using Propagator = Acts::Propagator<MultiStepper, Acts::Navigator>;
using DirectPropagator = Acts::Propagator<MultiStepper, Acts::DirectNavigator>;

using Fitter =
    Acts::Experimental::GaussianSumFitter<Propagator, BetheHeitlerApprox,
                                          Acts::VectorMultiTrajectory>;
using DirectFitter =
    Acts::Experimental::GaussianSumFitter<DirectPropagator, BetheHeitlerApprox,
                                          Acts::VectorMultiTrajectory>;
using TrackContainer =
    Acts::TrackContainer<Acts::VectorTrackContainer,
                         Acts::VectorMultiTrajectory, std::shared_ptr>;

struct GsfFitterFunctionImpl final : public ActsExamples::TrackFitterFunction {
  Fitter fitter;
  DirectFitter directFitter;

  Acts::GainMatrixUpdater updater;

  std::size_t maxComponents = 0;
  double weightCutoff = 0;
  bool abortOnError = false;
  bool disableAllMaterialHandling = false;

  GsfFitterFunctionImpl(Fitter&& f, DirectFitter&& df)
      : fitter(std::move(f)), directFitter(std::move(df)) {}

  template <typename calibrator_t>
  auto makeGsfOptions(const GeneralFitterOptions& options,
                      const calibrator_t& calibrator) const {
    Acts::Experimental::GsfExtensions<Acts::VectorMultiTrajectory> extensions;
    extensions.updater.connect<
        &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
        &updater);

    Acts::Experimental::GsfOptions<Acts::VectorMultiTrajectory> gsfOptions{
        options.geoContext,
        options.magFieldContext,
        options.calibrationContext,
        extensions,
        options.propOptions,
        &(*options.referenceSurface),
        maxComponents,
        weightCutoff,
        abortOnError,
        disableAllMaterialHandling};

    gsfOptions.extensions.calibrator.connect<&calibrator_t::calibrate>(
        &calibrator);

    return gsfOptions;
  }

  TrackFitterResult operator()(const std::vector<Acts::SourceLink>& sourceLinks,
                               const TrackParameters& initialParameters,
                               const GeneralFitterOptions& options,
                               const MeasurementCalibrator& calibrator,
                               TrackContainer& tracks) const override {
    const auto gsfOptions = makeGsfOptions(options, calibrator);

    using namespace Acts::Experimental::GsfConstants;
    if (not tracks.hasColumn(
            Acts::hashString(kFinalMultiComponentStateColumn))) {
      std::string key(kFinalMultiComponentStateColumn);
      tracks.template addColumn<FinalMultiComponentState>(key);
    }

    return fitter.fit(sourceLinks.begin(), sourceLinks.end(), initialParameters,
                      gsfOptions, tracks);
  }

  TrackFitterResult operator()(
      const std::vector<Acts::SourceLink>& sourceLinks,
      const TrackParameters& initialParameters,
      const GeneralFitterOptions& options,
      const RefittingCalibrator& calibrator,
      const std::vector<const Acts::Surface*>& surfaceSequence,
      TrackContainer& tracks) const override {
    const auto gsfOptions = makeGsfOptions(options, calibrator);

    using namespace Acts::Experimental::GsfConstants;
    if (not tracks.hasColumn(
            Acts::hashString(kFinalMultiComponentStateColumn))) {
      std::string key(kFinalMultiComponentStateColumn);
      tracks.template addColumn<FinalMultiComponentState>(key);
    }

    return directFitter.fit(sourceLinks.begin(), sourceLinks.end(),
                            initialParameters, gsfOptions, surfaceSequence,
                            tracks);
  }
};

}  // namespace

std::shared_ptr<TrackFitterFunction> ActsExamples::makeGsfFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    BetheHeitlerApprox betheHeitlerApprox, std::size_t maxComponents,
    double weightCutoff, Acts::FinalReductionMethod finalReductionMethod,
    bool abortOnError, bool disableAllMaterialHandling,
    const Acts::Logger& logger) {
  // Standard fitter
  MultiStepper stepper(magneticField, finalReductionMethod,
                       logger.cloneWithSuffix("Step"));
  Acts::Navigator::Config cfg{std::move(trackingGeometry)};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg, logger.cloneWithSuffix("Navigator"));
  Propagator propagator(std::move(stepper), std::move(navigator),
                        logger.cloneWithSuffix("Propagator"));
  Fitter trackFitter(std::move(propagator),
                     BetheHeitlerApprox(betheHeitlerApprox),
                     logger.cloneWithSuffix("GSF"));

  // Direct fitter
  MultiStepper directStepper(std::move(magneticField), finalReductionMethod,
                             logger.cloneWithSuffix("Step"));
  Acts::DirectNavigator directNavigator{
      logger.cloneWithSuffix("DirectNavigator")};
  DirectPropagator directPropagator(std::move(directStepper),
                                    std::move(directNavigator),
                                    logger.cloneWithSuffix("DirectPropagator"));
  DirectFitter directTrackFitter(std::move(directPropagator),
                                 BetheHeitlerApprox(betheHeitlerApprox),
                                 logger.cloneWithSuffix("DirectGSF"));

  // build the fitter functions. owns the fitter object.
  auto fitterFunction = std::make_shared<GsfFitterFunctionImpl>(
      std::move(trackFitter), std::move(directTrackFitter));
  fitterFunction->maxComponents = maxComponents;
  fitterFunction->weightCutoff = weightCutoff;
  fitterFunction->abortOnError = abortOnError;
  fitterFunction->disableAllMaterialHandling = disableAllMaterialHandling;

  return fitterFunction;
}
