// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief  Task for analysis of rho in UPCs using UD tables (from SG producer).
/// \author Edgardo Olivares, edgardo.olivares@alumno.buap.mx

#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
#include "TH1F.h"
#include "TH2F.h"

#include "random"
#include <numeric>
#include <string>
#include <vector>
#include <cmath>

using namespace o2;
using namespace o2::framework;

// Define UD tables
using UDtracks = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
using UDCollisions = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDZdcsReduced>;

namespace o2::aod
{
namespace twopi
{
// Declare columns
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);                        // Run number for event identification
DECLARE_SOA_COLUMN(M, m, double);                                         // Invariant mass of the system
DECLARE_SOA_COLUMN(Pt, pt, double);                                       // Transverse momentum of the system
DECLARE_SOA_COLUMN(Eta, eta, double);                                     // Pseudorapidity of the system
DECLARE_SOA_COLUMN(Phi, phi, double);                                     // Azimuthal angle of the system
DECLARE_SOA_COLUMN(PosX, posX, double);                                   // Vertex X position
DECLARE_SOA_COLUMN(PosY, posY, double);                                   // Vertex Y position
DECLARE_SOA_COLUMN(PosZ, posZ, double);                                   // Vertex Z position
DECLARE_SOA_COLUMN(NumContrib, numContrib, int32_t);                      // Number of primary vertex contributors

} // namespace twopi

// Define the output
DECLARE_SOA_TABLE(SYSTEMTREE, "AOD", "SystemTree",
                  twopi::RunNumber, twopi::M, twopi::Pt, twopi::Eta, twopi::Phi,
                  twopi::PosX, twopi::PosY, twopi::PosZ, twopi::NumContrib); 
} // namespace o2::aod

    struct upcTwoPionAnalysis {
  Produces<aod::SYSTEMTREE> systemTree;

  // System selection configuration
  Configurable<double> systemYCut{"systemYCut", 0.5, "Max Rapidity of rho prime"};
  //Configurable<double> systemPtCut{"systemPtCut", 0.1, "Min Pt of rho prime"};
  Configurable<double> systemMassMinCut{"systemMassMinCut", 0.0, "Min Mass of rho prime"};
  Configurable<double> systemMassMaxCut{"systemMassMaxCut", 1.2, "Max Mass of rho prime"};
  Configurable<double> etaCut{"etaCut", 0.9, "Track Pseudorapidity"};
  Configurable<double> numPVContrib{"numPVContrib", 2, "Number of PV contributors"};

  // AÑADIDO: Aux funct
  static double eta(double px, double py, double pz) {
    double p = std::sqrt(px*px + py*py + pz*pz);
    return 0.5 * std::log((p + pz) / (p - pz));
  }

  static double phi(double px, double py) {
    return std::atan2(py, px);
  }


  void process(UDCollisions::iterator const& collision, UDtracks const& tracks)
  {
    if (collision.numContrib() != numPVContrib)
    	return;
    
    std::vector<decltype(tracks.begin())> posPions;
    std::vector<decltype(tracks.begin())> negPions;

    // Loop over all tracks in the event
    for (const auto& track : tracks) {

      float trackEta = eta(track.px(), track.py(), track.pz());

      if (std::abs(trackEta) > etaCut) {
        continue;
      }

      if (std::abs(track.tpcNSigmaPi()) > 5.0) {
        continue;
      }

      // Track passed all selection criteria

      if (track.sign() > 0 && posPions.size() < 1) {
        posPions.push_back(track);
      } else if (track.sign() < 0 && negPions.size() < 1) {
        negPions.push_back(track);
      }

      if (posPions.size() == 1 && negPions.size() == 1)
        break;
    }

    if (posPions.size() != 1 || negPions.size() != 1) {
      return;
    }

    // Reconstruct the 2-pion system
    ROOT::Math::PxPyPzMVector twoPionSystem;

    for (const auto& track : {posPions[0], negPions[0]}) {
      ROOT::Math::PxPyPzMVector pionVec(
        track.px(), track.py(), track.pz(),
        o2::constants::physics::MassPionCharged);
      twoPionSystem += pionVec;
    }
    
   // Apply system-level kinematic cuts
    if (twoPionSystem.M() < systemMassMinCut || twoPionSystem.M() > systemMassMaxCut)
      return;
   // if (twoPionSystem.Pt() > systemPtCut)
     // return;
    if (std::abs(twoPionSystem.Rapidity()) > systemYCut)
      return;

    // Fill the output
    systemTree(
      collision.runNumber(),
      twoPionSystem.M(),
      twoPionSystem.Pt(),
      twoPionSystem.Rapidity(),
      twoPionSystem.Phi(),
      collision.posX(),
      collision.posY(),
      collision.posZ(),
      collision.numContrib());
   }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<upcTwoPionAnalysis>(cfgc, TaskName{"upc-two-pion-analysis"})};
}
