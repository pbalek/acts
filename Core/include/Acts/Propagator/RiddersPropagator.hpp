// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/PropagatorOptions.hpp"
#include "Acts/Propagator/PropagatorResult.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

template <typename propagator_t, typename actor_list_t = ActorList<>>
struct RiddersPropagatorOptions
    : public propagator_t::template Options<actor_list_t> {
  using base_type = propagator_t::template Options<actor_list_t>;

  using stepper_options_type = typename base_type::stepper_options_type;
  using navigator_options_type = typename base_type::navigator_options_type;
  using actor_list_type = actor_list_t;

  /// PropagatorOptions with context
  RiddersPropagatorOptions(const GeometryContext& gctx,
                           const MagneticFieldContext& mctx)
      : base_type(gctx, mctx) {}

  /// PropagatorOptions with context and plain options
  explicit RiddersPropagatorOptions(const PropagatorPlainOptions& pOptions)
      : base_type(pOptions) {}

  using base_type::operator PropagatorPlainOptions;

  /// @brief Expand the options with extended actors
  ///
  /// @tparam extended_actor_list_t Type of the new actor list
  ///
  /// @param extendedActorList The new actor list to be used (internally)
  template <typename extended_actor_list_t>
  RiddersPropagatorOptions<propagator_t, extended_actor_list_t> extend(
      extended_actor_list_t extendedActorList) const {
    RiddersPropagatorOptions<propagator_t, extended_actor_list_t> eoptions(
        base_type::geoContext, base_type::magFieldContext);

    static_cast<decltype(eoptions)::base_type&>(eoptions) =
        base_type::extend(std::move(extendedActorList));

    return eoptions;
  }

  using base_type::setPlainOptions;

  /// Initial scale for the deviation of the individual bound track parameters
  BoundVector deviationScale = {1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4};

  /// Different factors applied to the initial scale to create the
  /// deviations of the individual bound track parameters. The resulting
  /// function value deviations are then fitted to a line to determine the
  /// first order derivatives of the final parameters wrt. the initial
  /// parameters.
  std::vector<double> deviationFactors = {-2, -1, 1, 2};
};

/// @brief This class performs the Ridders algorithm to estimate the propagation
/// of the covariance to a certain point in space.
///
/// The algorithm is based on the small deviations of the start parameters based
/// on their uncertainty at the beginning of the propgation. This deviation is
/// represented here by a vector of relative deviations of these parameters and
/// fix for all parameters. So, a common choice has to be found that is able to
/// actually fit into the order of magnitude of the uncertainty of each
/// parameter. Using these deviations, the propagation is repeated multiple
/// times and the final covariance matrix at a given target surface is
/// afterwards evaluated by first order derivatives of the final state
/// parameters wrt. the initial parameters. Therefore this evaluation represents
/// a first order approximation of the transport jacobian. Since performing
/// multiple propagations and a numerical evaluation of the covariance requires
/// more time than a single propagation towards a target + a common propagation
/// of the covariance, this class just serves to verify the results of the
/// latter classes.
template <typename propagator_t>
class RiddersPropagator {
  ///
  /// @note The result_type_helper struct and the actor_list_t_result_t are
  /// here to allow a look'n'feel of this class like the Propagator itself
  ///

  /// @brief Helper struct determining the result's type
  ///
  /// @tparam parameters_t Type of final track parameters
  /// @tparam actor_list_t    List of propagation action types
  ///
  /// This helper struct provides type definitions to extract the correct
  /// propagation result type from a given TrackParameter type and an
  /// ActorList.
  ///
  template <typename parameters_t, typename actor_list_t>
  struct result_type_helper {
    /// @brief Propagation result type for an arbitrary list of additional
    ///        propagation results
    ///
    /// @tparam args Parameter pack specifying additional propagation results
    ///
    template <typename... args>
    using this_result_type = PropagatorResult<parameters_t, args...>;

    /// @brief Propagation result type derived from a given action list
    using type = typename actor_list_t::template result_type<this_result_type>;
  };

  /// @brief Short-hand type definition for propagation result derived from
  ///        an action list
  ///
  /// @tparam parameters_t Type of the final track parameters
  /// @tparam actor_list_t List of propagation action types
  ///
  template <typename parameters_t, typename actor_list_t>
  using actor_list_t_result_t =
      typename result_type_helper<parameters_t, actor_list_t>::type;

 public:
  /// Type of the stepper in use for public scope
  using Stepper = typename propagator_t::Stepper;

  /// Type of the navigator in use for public scope
  using Navigator = typename propagator_t::Navigator;

  /// Type of state object used by the propagation implementation
  using StepperState = typename Stepper::State;

  /// Typedef the navigator state
  using NavigatorState = typename Navigator::State;

  template <typename propagator_options_t, typename... extension_state_t>
  using State = typename propagator_t::template State<
      propagator_options_t, StepperState, NavigatorState, extension_state_t...>;

  using StepperOptions = typename Stepper::Options;

  using NavigatorOptions = typename Navigator::Options;

  template <typename actor_list_t = ActorList<>>
  using Options = RiddersPropagatorOptions<propagator_t, actor_list_t>;

  /// @brief Constructor using a propagator
  ///
  /// @param [in] propagator The propagator to use
  explicit RiddersPropagator(propagator_t propagator)
      : m_propagator(std::move(propagator)) {}

  /// @brief Propagation method targeting curvilinear parameters
  ///
  /// @tparam parameters_t Type of the start parameters
  /// @tparam propagator_options_t Type of the propagator options
  ///
  /// @param [in] start Start parameters
  /// @param [in] options Options of the propagations
  ///
  /// @return Result of the propagation
  template <typename parameters_t, typename propagator_options_t>
  Result<actor_list_t_result_t<BoundTrackParameters,
                               typename propagator_options_t::actor_list_type>>
  propagate(const parameters_t& start,
            const propagator_options_t& options) const;

  /// @brief Propagation method targeting bound parameters
  ///
  /// @tparam parameters_t Type of the start parameters
  /// @tparam propagator_options_t Type of the propagator options
  ///
  /// @param [in] start Start parameters
  /// @param [in] target The target surface
  /// @param [in] options Options of the propagations
  ///
  /// @return Result of the propagation
  /// @note If the target surface is a disc, the resulting covariance may be
  /// inconsistent. In this case a zero matrix is returned.
  template <typename parameters_t, typename propagator_options_t>
  Result<actor_list_t_result_t<BoundTrackParameters,
                               typename propagator_options_t::actor_list_type>>
  propagate(const parameters_t& start, const Surface& target,
            const propagator_options_t& options) const;

 private:
  /// Does the actual ridders propagation by wiggling the parameters and
  /// propagating again. This function is called from the different
  /// propagation overloads in order to deduplicate code.
  ///
  /// @param [in] options Options of the propagations
  /// @param [in] start Start parameters
  /// @param [in] nominalResult The result of the nominal propagation
  template <typename parameters_t, typename propagator_options_t>
  BoundMatrix wiggleAndCalculateJacobian(
      const parameters_t& start, const propagator_options_t& options,
      const actor_list_t_result_t<
          BoundTrackParameters, typename propagator_options_t::actor_list_type>&
          nominalResult) const;

  /// @brief This function wiggles one dimension of the starting parameters,
  /// performs the propagation to a surface and collects for each change of the
  /// start parameters the slope
  ///
  /// @tparam options_t PropagatorOptions object
  /// @tparam parameters_t Type of the parameters to start the propagation with
  ///
  /// @param [in] options Options do define how to wiggle
  /// @param [in] start Start parameters which will be modified
  /// @param [in] param Index to get the parameter that will be modified
  /// @param [in] target Target surface
  /// @param [in] nominal Nominal end parameters
  /// @param [in] deviations Vector of deviations
  ///
  /// @return Vector containing each slope
  template <typename propagator_options_t, typename parameters_t>
  std::vector<BoundVector> wiggleParameter(
      const propagator_options_t& options, const parameters_t& start,
      unsigned int param, const Surface& target, const BoundVector& nominal,
      const std::vector<double>& deviations) const;

  propagator_t m_propagator;
};

}  // namespace Acts

#include "Acts/Propagator/RiddersPropagator.ipp"
