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
/// \brief  Task for analysis of rho prime in UPCs using UD tables (from SG producer).
/// \author Cesar Ramirez, cesar.ramirez@cern.ch

#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/PIDResponse.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
#include "TH1F.h"
#include "TH2F.h"

#include "random"
#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace ROOT::Math;

// Define UD tables
using UDtracks = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
using UDCollisions = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDZdcsReduced>;
using UDMcCollisions = aod::UDMcCollisions;
using UDMcParticles = aod::UDMcParticles;

namespace o2::aod
{
namespace fourpi
{
// Columns for real data
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);                        // Run number for event identification
DECLARE_SOA_COLUMN(M, m, double);                                         // Invariant mass of the system
DECLARE_SOA_COLUMN(Pt, pt, double);                                       // Transverse momentum of the system
DECLARE_SOA_COLUMN(Eta, eta, double);                                     // Pseudorapidity of the system
DECLARE_SOA_COLUMN(Phi, phi, double);                                     // Azimuthal angle of the system
DECLARE_SOA_COLUMN(PosX, posX, double);                                   // Vertex X position
DECLARE_SOA_COLUMN(PosY, posY, double);                                   // Vertex Y position
DECLARE_SOA_COLUMN(PosZ, posZ, double);                                   // Vertex Z position
DECLARE_SOA_COLUMN(TotalCharge, totalCharge, int);                        // Total charge of selected tracks
DECLARE_SOA_COLUMN(TotalFT0AmplitudeA, totalFT0AmplitudeA, float);        // FT0A amplitude
DECLARE_SOA_COLUMN(TotalFT0AmplitudeC, totalFT0AmplitudeC, float);        // FT0C amplitude
DECLARE_SOA_COLUMN(TotalFV0AmplitudeA, totalFV0AmplitudeA, float);        // FV0A amplitude
DECLARE_SOA_COLUMN(NumContrib, numContrib, int32_t);                      // Number of primary vertex contributors
DECLARE_SOA_COLUMN(Sign, sign, std::vector<int>);                         // Track charges
DECLARE_SOA_COLUMN(TrackPt, trackPt, std::vector<float>);                 // Track pT values
DECLARE_SOA_COLUMN(TrackEta, trackEta, std::vector<float>);               // Track eta values
DECLARE_SOA_COLUMN(TrackPhi, trackPhi, std::vector<float>);               // Track phi values
DECLARE_SOA_COLUMN(TPCNSigmaEl, tpcNSigmaEl, std::vector<float>);         // TPC nσ for electrons
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, std::vector<float>);         // TPC nσ for pions
DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcNSigmaKa, std::vector<float>);         // TPC nσ for kaons
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, std::vector<float>);         // TPC nσ for protons
DECLARE_SOA_COLUMN(TrackID, trackID, std::vector<int>);                   // Track identifiers
DECLARE_SOA_COLUMN(IsReconstructedWithUPC, isReconstructedWithUPC, bool); // UPC mode reconstruction flag
DECLARE_SOA_COLUMN(TimeZNA, timeZNA, float);                              // ZNA timing
DECLARE_SOA_COLUMN(TimeZNC, timeZNC, float);                              // ZNC timing
DECLARE_SOA_COLUMN(EnergyCommonZNA, energyCommonZNA, float);              // ZNA energy
DECLARE_SOA_COLUMN(EnergyCommonZNC, energyCommonZNC, float);              // ZNC energy
DECLARE_SOA_COLUMN(IsChargeZero, isChargeZero, bool);                     // Neutral system flag
DECLARE_SOA_COLUMN(OccupancyInTime, occupancyInTime, int);                // Occupancy in time
DECLARE_SOA_COLUMN(HadronicRate, hadronicRate, double);                   // Hadronic interaction rate
} // namespace fourpi

// Table for real data
DECLARE_SOA_TABLE(SYSTEMTREE, "AOD", "SystemTree",
                  fourpi::RunNumber, fourpi::M, fourpi::Pt, fourpi::Eta, fourpi::Phi,
                  fourpi::PosX, fourpi::PosY, fourpi::PosZ, fourpi::TotalCharge,
                  fourpi::TotalFT0AmplitudeA, fourpi::TotalFT0AmplitudeC, fourpi::TotalFV0AmplitudeA,
                  fourpi::NumContrib,
                  fourpi::Sign, fourpi::TrackPt, fourpi::TrackEta, fourpi::TrackPhi,
                  fourpi::TPCNSigmaEl, fourpi::TPCNSigmaPi, fourpi::TPCNSigmaKa, fourpi::TPCNSigmaPr,
                  fourpi::TrackID, fourpi::IsReconstructedWithUPC,
                  fourpi::TimeZNA, fourpi::TimeZNC, fourpi::EnergyCommonZNA, fourpi::EnergyCommonZNC,
                  fourpi::IsChargeZero, fourpi::OccupancyInTime, fourpi::HadronicRate);

namespace mcfourpi
{
// Columns for MC
DECLARE_SOA_COLUMN(RunNumberMC, runNumberMC, int32_t);
DECLARE_SOA_COLUMN(MCM, mcm, float);                  // True mass
DECLARE_SOA_COLUMN(MCPt, mcpt, float);                // True pT
DECLARE_SOA_COLUMN(MCEta, mceta, float);              // True eta
DECLARE_SOA_COLUMN(MCPhi, mcphi, float);              // True phi
DECLARE_SOA_COLUMN(MCPosX, mcposX, float);            // True vertex X
DECLARE_SOA_COLUMN(MCPosY, mcposY, float);            // True vertex Y
DECLARE_SOA_COLUMN(MCPosZ, mcposZ, float);            // True vertex Z
DECLARE_SOA_COLUMN(MCTrackSign, mctrackSign, int[4]); // Pion charges
DECLARE_SOA_COLUMN(MCTrackPt, mctrackPt, float[4]);   // Pion pTs
DECLARE_SOA_COLUMN(MCTrackEta, mctrackEta, float[4]); // Pion etas
DECLARE_SOA_COLUMN(MCTrackPhi, mctrackPhi, float[4]); // Pion phis
DECLARE_SOA_COLUMN(MCTrackPdg, mctrackPdg, int[4]);   // PDG codes
DECLARE_SOA_COLUMN(MCIsPhysicalPrimary0, mcIsPhysicalPrimary0, bool);
DECLARE_SOA_COLUMN(MCIsPhysicalPrimary1, mcIsPhysicalPrimary1, bool);
DECLARE_SOA_COLUMN(MCIsPhysicalPrimary2, mcIsPhysicalPrimary2, bool);
DECLARE_SOA_COLUMN(MCIsPhysicalPrimary3, mcIsPhysicalPrimary3, bool);
} // namespace mcfourpi

// Table for MC
DECLARE_SOA_TABLE(MCFourPiTree, "AOD", "MCFOURPITREE",
                  mcfourpi::RunNumberMC,
                  mcfourpi::MCM, mcfourpi::MCPt, mcfourpi::MCEta, mcfourpi::MCPhi,
                  mcfourpi::MCPosX, mcfourpi::MCPosY, mcfourpi::MCPosZ,
                  mcfourpi::MCTrackSign, mcfourpi::MCTrackPt,
                  mcfourpi::MCTrackEta, mcfourpi::MCTrackPhi,
                  mcfourpi::MCTrackPdg,
                  mcfourpi::MCIsPhysicalPrimary0, mcfourpi::MCIsPhysicalPrimary1,
                  mcfourpi::MCIsPhysicalPrimary2, mcfourpi::MCIsPhysicalPrimary3);
} // namespace o2::aod

struct upcRhoPrimeAnalysis {
  Produces<aod::SYSTEMTREE> systemTree;
  Produces<aod::MCFourPiTree> mcFourPiTree;

  SGSelector sgSelector;

  // System selection configuration
  Configurable<double> systemYCut{"systemYCut", 0.5, "Max Rapidity of rho prime"};
  Configurable<double> systemPtCut{"systemPtCut", 0.1, "Min Pt of rho prime"};
  Configurable<double> systemMassMinCut{"systemMassMinCut", 0.8, "Min Mass of rho prime"};
  Configurable<double> systemMassMaxCut{"systemMassMaxCut", 2.2, "Max Mass of rho prime"};
  Configurable<double> etaCut{"etaCut", 0.9, "Track Pseudorapidity"};

  // Event selection configuration
  Configurable<float> vZCut{"vZCut", 10.0, "Cut on vertex Z position"};
  Configurable<int> numPVContrib{"numPVContrib", 4, "Number of PV contributors"};
  Configurable<float> fv0Cut{"fv0Cut", 50.0, "FV0 amplitude cut"};
  Configurable<float> ft0aCut{"ft0aCut", 50.0, "FT0A amplitude cut"};
  Configurable<float> ft0cCut{"ft0cCut", 50.0, "FT0C amplitude cut"};
  Configurable<float> zdcCut{"zdcCut", 0.0, "ZDC energy cut"};
  Configurable<bool> sbpCut{"sbpCut", true, "SBP cut"};
  Configurable<bool> itsROFbCut{"itsROFbCut", true, "ITS ROFb cut"};
  Configurable<bool> vtxITSTPCcut{"vtxITSTPCcut", true, "Vertex ITS-TPC cut"};
  Configurable<bool> tfbCut{"tfbCut", true, "TFB cut"};
  Configurable<bool> specifyGapSide{"specifyGapSide", true, "specify gap side for SG/DG produced data"};
  Configurable<int> gapSide{"gapSide", 2, "gap side for SG produced data"};

  // Track selection configuration
  Configurable<bool> useOnlyPVtracks{"useOnlyPVtracks", true, "Use only PV tracks"};
  Configurable<float> tpcChi2NClsCut{"tpcChi2NClsCut", 5.0, "TPC chi2/N clusters cut"};
  Configurable<float> itsChi2NClsCut{"itsChi2NClsCut", 36.0, "ITS chi2/N clusters cut"};
  Configurable<float> nSigmaTPCcut{"nSigmaTPCcut", 5.0, "TPC nSigma cut"};
  Configurable<float> dcaXYcut{"dcaXYcut", 0, "dcaXY cut"};
  Configurable<float> dcaZcut{"dcaZcut", 2, "dcaZ cut"};
  Configurable<int> minTPCFindableClusters{"minTPCFindableClusters", 70, "Minimum number of findable TPC clusters"};

  Configurable<bool> cfgProcessMC{"cfgProcessMC", true, "Process MC data"};

  // Define histogram registry for real data analysis
  HistogramRegistry registry{
    "registry",
    {// Event flow histograms
     {"Events/Flow", "Event flow;Cut;Counts", {HistType::kTH1D, {{9, 0, 9}}}},
     {"Events/VertexZ", "Vertex Z;z (cm);Counts", {HistType::kTH1F, {{200, -20, 20}}}},
     {"Events/NumContrib", "Number of contributors;N_{contrib};Counts", {HistType::kTH1F, {{100, 0, 100}}}},
     {"Events/FV0Amplitude", "FV0 amplitude;Amplitude;Counts", {HistType::kTH1F, {{200, 0, 200}}}},
     {"Events/FT0AmplitudeA", "FT0A amplitude;Amplitude;Counts", {HistType::kTH1F, {{200, 0, 200}}}},
     {"Events/FT0AmplitudeC", "FT0C amplitude;Amplitude;Counts", {HistType::kTH1F, {{200, 0, 200}}}},
     {"Events/ZDCEnergy", "ZDC energy;Energy (TeV);Counts", {HistType::kTH1F, {{200, 0, 2}}}},

     // Track quality histograms
     {"Tracks/Pt", "Track p_{T};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{200, 0, 2}}}},
     {"Tracks/Eta", "Track #eta;#eta;Counts", {HistType::kTH1F, {{200, -2, 2}}}},
     {"Tracks/TPCNSigmaPi", "TPC n#sigma for #pi;n#sigma;Counts", {HistType::kTH1F, {{200, -10, 10}}}},
     {"Tracks/TPCChi2NCl", "TPC #chi^{2}/N_{cls};#chi^{2}/N_{cls};Counts", {HistType::kTH1F, {{200, 0, 20}}}},
     {"Tracks/ITSChi2NCl", "ITS #chi^{2}/N_{cls};#chi^{2}/N_{cls};Counts", {HistType::kTH1F, {{200, 0, 50}}}},
     {"Tracks/RejectionReasons", "Track rejection reasons;Reason;Counts", {HistType::kTH1F, {{12, 0, 12}}}},
     {"Tracks/DCASpectrum", "Track DCA spectrum;DCA (cm);Counts", {HistType::kTH1F, {{100, 0, 5}}}},
     {"Tracks/ChargeDistribution", "Track charge distribution;Charge;Counts", {HistType::kTH1F, {{3, -1.5, 1.5}}}},
     {"Tracks/TPCClusters", "TPC clusters findable;N_{clusters};Counts", {HistType::kTH1F, {{100, 0, 200}}}},

     // System kinematics histograms
     {"System/hM", ";m (GeV/#it{c}^{2});counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"System/hPt", ";p_{T} (GeV/#it{c});counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"System/hEta", ";#eta;counts", {HistType::kTH1F, {{180, -0.9, 0.9}}}},
     {"System/hPhi", ";#phi;counts", {HistType::kTH1F, {{180, 0.0, 6.28}}}},
     {"System/hY", ";y;counts", {HistType::kTH1F, {{180, -0.9, 0.9}}}},

     // Comparison histograms
     {"Cuts/MBefore", "Mass before cuts;m (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0, 10}}}},
     {"Cuts/MAfter", "Mass after cuts;m (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0, 10}}}},
     {"Cuts/PtBefore", "p_{T} before cuts;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0, 1}}}},
     {"Cuts/PtAfter", "p_{T} after cuts;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0, 10}}}}}};

  // Define histogram registry for MC analysis
  HistogramRegistry mcRegistry{
    "mcRegistry",
    {// Event-level MC histograms
     {"MC/Events/hAllEvents", "All MC Events", {HistType::kTH1F, {{1, 0, 1}}}},
     {"MC/Events/hAccepted", "Accepted MC Events", {HistType::kTH1F, {{1, 0, 1}}}},
     {"MC/Events/hVertexZ", "MC Vertex Z;z (cm);Counts", {HistType::kTH1F, {{400, -20.0, 20.0}}}},
     {"MC/Events/hNPions", "Number of primary pions;N_{#pi};Counts", {HistType::kTH1F, {{10, -0.5, 9.5}}}},

     // Track-level MC histograms
     {"MC/Tracks/hPt", "MC Track p_{T};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{200, 0.0, 2.0}}}},
     {"MC/Tracks/hEta", "MC Track #eta;#eta;Counts", {HistType::kTH1F, {{200, -2.0, 2.0}}}},
     {"MC/Tracks/hPhi", "MC Track #phi;#phi;Counts", {HistType::kTH1F, {{200, 0.0, 6.28}}}},
     {"MC/Tracks/hPdgCode", "PDG codes;PDG code;Counts", {HistType::kTH1F, {{6001, -3000.5, 3000.5}}}},

     // System-level MC histograms
     {"MC/System/hM", "MC Invariant Mass;m (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"MC/System/hPt", "MC p_{T};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"MC/System/hY", "MC Rapidity;y;Counts", {HistType::kTH1F, {{180, -0.9, 0.9}}}},
     {"MC/System/hMvsPt", "MC Mass vs p_{T};m (GeV/c^{2});p_{T} (GeV/c)", {HistType::kTH2F, {{1000, 0.0, 10.0}, {1000, 0.0, 10.0}}}},
     {"MC/System/hMvsY", "MC Mass vs Rapidity;m (GeV/c^{2});y", {HistType::kTH2F, {{1000, 0.0, 10.0}, {200, -2.0, 2.0}}}},

     // Rho prime specific histograms
     {"MC/RhoPrime/hFound", "Rho Prime Found;Found;Counts", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
     {"MC/RhoPrime/hMass", "Rho Prime Mass;m (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"MC/RhoPrime/hPt", "Rho Prime p_{T};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"MC/RhoPrime/hDecayPions", "Rho Prime Decay Pions;N_{#pi};Counts", {HistType::kTH1F, {{5, -0.5, 4.5}}}},

     // Control histograms for MC analysis
     {"MC/Control/hNDaughters", "Number of Daughters;N_{daughters};Counts", {HistType::kTH1F, {{10, 0, 10}}}},
     {"MC/Control/hNPiPlus", "Number of #pi^{+};N_{#pi^{+}};Counts", {HistType::kTH1F, {{5, -0.5, 4.5}}}},
     {"MC/Control/hNPiMinus", "Number of #pi^{-};N_{#pi^{-}};Counts", {HistType::kTH1F, {{5, -0.5, 4.5}}}},
     {"MC/Control/hMAll", "MC Mass All;m (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"MC/Control/hPtAll", "MC p_{T} All;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"MC/Control/hYAll", "MC Rapidity All;y;Counts", {HistType::kTH1F, {{200, -2.0, 2.0}}}},
     {"MC/Control/hMvsPtAll", "MC Mass vs p_{T} All;m (GeV/c^{2});p_{T} (GeV/c)", {HistType::kTH2F, {{500, 0.0, 5.0}, {500, 0.0, 5.0}}}},

     // Cut rejection histograms for MC
     {"MC/Cuts/hMassRejected", "Rejected by Mass;m (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"MC/Cuts/hPtRejected", "Rejected by p_{T};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"MC/Cuts/hYRejected", "Rejected by Rapidity;y;Counts", {HistType::kTH1F, {{200, -2.0, 2.0}}}}}};

  // Helper functions for kinematic calculations
  static float pt(float px, float py) { return std::sqrt(px * px + py * py); }

  static float eta(float px, float py, float pz)
  {
    if (std::abs(pz) > 1e-10) {
      float p = std::sqrt(px * px + py * py + pz * pz);
      return 0.5f * std::log((p + pz) / (p - pz));
    }
    return 0.0f;
  }

  static float phi(float px, float py)
  {
    if (std::abs(px) > 1e-10 || std::abs(py) > 1e-10) {
      return std::atan2(py, px);
    }
    return 0.0f;
  }

  void init(InitContext&)
  {
    // Configure event flow histogram labels
    auto hFlow = registry.get<TH1>(HIST("Events/Flow"));
    hFlow->GetXaxis()->SetBinLabel(1, "All events");
    hFlow->GetXaxis()->SetBinLabel(2, "ITS-TPC cut");
    hFlow->GetXaxis()->SetBinLabel(3, "SBP cut");
    hFlow->GetXaxis()->SetBinLabel(4, "ITS ROFb cut");
    hFlow->GetXaxis()->SetBinLabel(5, "TFB cut");
    hFlow->GetXaxis()->SetBinLabel(6, "Gap Side cut");
    hFlow->GetXaxis()->SetBinLabel(7, "PV contrib cut");
    hFlow->GetXaxis()->SetBinLabel(8, "Z vtx cut");
    hFlow->GetXaxis()->SetBinLabel(9, "4 tracks cut");

    // Configure track rejection reasons histogram labels
    auto hReject = registry.get<TH1>(HIST("Tracks/RejectionReasons"));
    hReject->GetXaxis()->SetBinLabel(1, "All Tracks");
    hReject->GetXaxis()->SetBinLabel(2, "PV Contributor");
    hReject->GetXaxis()->SetBinLabel(3, "Has ITS+TPC");
    hReject->GetXaxis()->SetBinLabel(4, "pT > 0.1 GeV/c");
    hReject->GetXaxis()->SetBinLabel(5, "TPC chi2/cluster");
    hReject->GetXaxis()->SetBinLabel(6, "ITS chi2/cluster");
    hReject->GetXaxis()->SetBinLabel(7, "TPC clusters findable");
    hReject->GetXaxis()->SetBinLabel(8, "TPC nSigmaPi");
    hReject->GetXaxis()->SetBinLabel(9, "Eta acceptance");
    hReject->GetXaxis()->SetBinLabel(10, "DCAz cut");
    hReject->GetXaxis()->SetBinLabel(11, "DCAxy cut");
    hReject->GetXaxis()->SetBinLabel(12, "Accepted Tracks");

    // Configure labels for MC histograms
    auto hPdg = mcRegistry.get<TH1>(HIST("MC/Tracks/hPdgCode"));
    hPdg->GetXaxis()->SetBinLabel(hPdg->GetXaxis()->FindBin(211), "#pi^{+}");
    hPdg->GetXaxis()->SetBinLabel(hPdg->GetXaxis()->FindBin(-211), "#pi^{-}");
    hPdg->GetXaxis()->SetBinLabel(hPdg->GetXaxis()->FindBin(30113), "rho'");

    auto hRhoFound = mcRegistry.get<TH1>(HIST("MC/RhoPrime/hFound"));
    hRhoFound->GetXaxis()->SetBinLabel(1, "Not Found");
    hRhoFound->GetXaxis()->SetBinLabel(2, "Found");
  }

  // Process for real data
  void processRealData(UDCollisions::iterator const& collision, UDtracks const& tracks)
  {
    // Count all processed events
    registry.fill(HIST("Events/Flow"), 0);

    // Fill basic event diagnostics
    registry.fill(HIST("Events/VertexZ"), collision.posZ());
    registry.fill(HIST("Events/NumContrib"), collision.numContrib());
    registry.fill(HIST("Events/FV0Amplitude"), collision.totalFV0AmplitudeA());
    registry.fill(HIST("Events/FT0AmplitudeA"), collision.totalFT0AmplitudeA());
    registry.fill(HIST("Events/FT0AmplitudeC"), collision.totalFT0AmplitudeC());
    registry.fill(HIST("Events/ZDCEnergy"), collision.energyCommonZNA());
    registry.fill(HIST("Events/ZDCEnergy"), collision.energyCommonZNC());

    // Apply event selection cuts in sequence
    if (collision.vtxITSTPC() != vtxITSTPCcut)
      return;
    registry.fill(HIST("Events/Flow"), 1);

    if (collision.sbp() != sbpCut)
      return;
    registry.fill(HIST("Events/Flow"), 2);

    if (collision.itsROFb() != itsROFbCut)
      return;
    registry.fill(HIST("Events/Flow"), 3);

    if (collision.tfb() != tfbCut)
      return;
    registry.fill(HIST("Events/Flow"), 4);

    if (specifyGapSide && collision.gapSide() != gapSide)
      return;
    if (collision.totalFV0AmplitudeA() > fv0Cut)
      return;
    if (collision.totalFT0AmplitudeA() > ft0aCut)
      return;
    if (collision.totalFT0AmplitudeC() > ft0cCut)
      return;
    if (collision.energyCommonZNA() > zdcCut || collision.energyCommonZNC() > zdcCut)
      return;
    registry.fill(HIST("Events/Flow"), 5);

    if (collision.numContrib() != numPVContrib)
      return;
    registry.fill(HIST("Events/Flow"), 6);

    if (std::abs(collision.posZ()) > vZCut)
      return;
    registry.fill(HIST("Events/Flow"), 7);

    // Select positive and negative pions
    std::vector<decltype(tracks.begin())> posPions;
    std::vector<decltype(tracks.begin())> negPions;
    posPions.reserve(2);
    negPions.reserve(2);

    // Loop over all tracks in the event
    for (const auto& track : tracks) {
      registry.fill(HIST("Tracks/RejectionReasons"), 0); // Count all tracks

      // Apply track selection criteria in sequence
      if (useOnlyPVtracks && !track.isPVContributor()) {
        registry.fill(HIST("Tracks/RejectionReasons"), 1);
        continue;
      }

      if (!track.hasITS() || !track.hasTPC()) {
        registry.fill(HIST("Tracks/RejectionReasons"), 2);
        continue;
      }

      // Fill track quality histograms
      registry.fill(HIST("Tracks/Pt"), track.pt());
      registry.fill(HIST("Tracks/Eta"), eta(track.px(), track.py(), track.pz()));
      registry.fill(HIST("Tracks/TPCNSigmaPi"), track.tpcNSigmaPi());
      registry.fill(HIST("Tracks/TPCChi2NCl"), track.tpcChi2NCl());
      registry.fill(HIST("Tracks/ITSChi2NCl"), track.itsChi2NCl());
      registry.fill(HIST("Tracks/DCASpectrum"), std::hypot(track.dcaXY(), track.dcaZ()));
      registry.fill(HIST("Tracks/ChargeDistribution"), track.sign());
      registry.fill(HIST("Tracks/TPCClusters"), track.tpcNClsFindable());

      if (track.pt() <= 0.1f) {
        registry.fill(HIST("Tracks/RejectionReasons"), 3);
        continue;
      }

      if (track.tpcChi2NCl() > tpcChi2NClsCut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 4);
        continue;
      }
      if (track.itsChi2NCl() > itsChi2NClsCut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 5);
        continue;
      }

      if (track.tpcNClsFindable() < minTPCFindableClusters) {
        registry.fill(HIST("Tracks/RejectionReasons"), 6);
        continue;
      }

      if (std::abs(track.tpcNSigmaPi()) > nSigmaTPCcut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 7);
        continue;
      }

      float trackEta = eta(track.px(), track.py(), track.pz());
      if (std::abs(trackEta) > etaCut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 8);
        continue;
      }

      if (std::abs(track.dcaZ()) > dcaZcut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 9);
        continue;
      }

      float maxDCAxy = 0.0105 + 0.035 / std::pow(track.pt(), 1.1);
      if (dcaXYcut == 0 && (std::fabs(track.dcaXY())) > maxDCAxy) {
        registry.fill(HIST("Tracks/RejectionReasons"), 10);
        continue;
      } else if (dcaXYcut == 0 && (std::fabs(track.dcaXY()) > maxDCAxy)) {
        registry.fill(HIST("Tracks/RejectionReasons"), 10);
        continue;
      }

      // Track passed all selection criteria
      registry.fill(HIST("Tracks/RejectionReasons"), 11);

      // Add to positive or negative pion lists
      if (track.sign() > 0 && posPions.size() < 2) {
        posPions.push_back(track);
      } else if (track.sign() < 0 && negPions.size() < 2) {
        negPions.push_back(track);
      }

      // Stop if we have enough tracks
      if (posPions.size() == 2 && negPions.size() == 2)
        break;
    }

    // Require exactly 2 positive and 2 negative pions
    if (posPions.size() != 2 || negPions.size() != 2) {
      return;
    }
    registry.fill(HIST("Events/Flow"), 8);

    // Combine selected tracks
    std::vector<decltype(tracks.begin())> selectedTracks;
    selectedTracks.insert(selectedTracks.end(), posPions.begin(), posPions.end());
    selectedTracks.insert(selectedTracks.end(), negPions.begin(), negPions.end());

    // Reconstruct the 4-pion system
    PxPyPzMVector fourPionSystem;
    std::vector<PxPyPzMVector> pionFourVectors;

    for (const auto& track : selectedTracks) {
      PxPyPzMVector pionVec(
        track.px(), track.py(), track.pz(),
        o2::constants::physics::MassPionCharged);
      fourPionSystem += pionVec;
      pionFourVectors.push_back(pionVec);
    }

    // Fill pre-cut system histograms
    registry.fill(HIST("Cuts/MBefore"), fourPionSystem.M());
    registry.fill(HIST("Cuts/PtBefore"), fourPionSystem.Pt());

    // Apply system-level kinematic cuts
    if (fourPionSystem.M() < systemMassMinCut || fourPionSystem.M() > systemMassMaxCut)
      return;
    if (fourPionSystem.Pt() > systemPtCut)
      return;
    if (std::abs(fourPionSystem.Rapidity()) > systemYCut)
      return;

    // Fill post-cut system histograms
    registry.fill(HIST("Cuts/MAfter"), fourPionSystem.M());
    registry.fill(HIST("Cuts/PtAfter"), fourPionSystem.Pt());
    registry.fill(HIST("System/hM"), fourPionSystem.M());
    registry.fill(HIST("System/hPt"), fourPionSystem.Pt());
    registry.fill(HIST("System/hEta"), fourPionSystem.Eta());
    registry.fill(HIST("System/hPhi"), fourPionSystem.Phi() + o2::constants::math::PI);
    registry.fill(HIST("System/hY"), fourPionSystem.Rapidity());

    // Prepare track information for output tree
    std::vector<float> trackPts, trackEtas, trackPhis;
    std::vector<int> trackSigns, trackIDs;
    std::vector<float> tpcNSigmasEl, tpcNSigmasPi, tpcNSigmasKa, tpcNSigmasPr;

    for (size_t i = 0; i < selectedTracks.size(); i++) {
      const auto& track = selectedTracks[i];
      trackPts.push_back(track.pt());
      trackEtas.push_back(eta(track.px(), track.py(), track.pz()));
      trackPhis.push_back(phi(track.px(), track.py()));
      trackSigns.push_back(track.sign());
      tpcNSigmasEl.push_back(track.tpcNSigmaEl());
      tpcNSigmasPi.push_back(track.tpcNSigmaPi());
      tpcNSigmasKa.push_back(track.tpcNSigmaKa());
      tpcNSigmasPr.push_back(track.tpcNSigmaPr());
      trackIDs.push_back(i);
    }

    bool isReconstructedWithUPC = (collision.flags() == 1);

    // Fill the output tree
    systemTree(
      collision.runNumber(),
      fourPionSystem.M(),
      fourPionSystem.Pt(),
      fourPionSystem.Rapidity(),
      fourPionSystem.Phi(),
      collision.posX(),
      collision.posY(),
      collision.posZ(),
      0, // Total charge = 0
      collision.totalFT0AmplitudeA(),
      collision.totalFT0AmplitudeC(),
      collision.totalFV0AmplitudeA(),
      collision.numContrib(),
      trackSigns,
      trackPts,
      trackEtas,
      trackPhis,
      tpcNSigmasEl,
      tpcNSigmasPi,
      tpcNSigmasKa,
      tpcNSigmasPr,
      trackIDs,
      isReconstructedWithUPC,
      collision.timeZNA(),
      collision.timeZNC(),
      collision.energyCommonZNA(),
      collision.energyCommonZNC(),
      true, // Always charge zero
      collision.occupancyInTime(),
      collision.hadronicRate());
  }

  // Process for MC generated data
void processMCgen(UDMcCollisions::iterator const& mcCollision, UDMcParticles const& mcParticles)
{
    // OBTENER RUN NUMBER DESDE GlobalBC - SIN PARÉNTESIS
    uint64_t globalBC = mcCollision.GlobalBC;  // ← SIN PARÉNTESIS
    int runNumber = (globalBC >> 32) & 0xFFFFFFFF;
    int runIndex = runNumber;
    // Constants for the analysis
    const int numFourPionTracks = 4;
    const int numPiPlus = 2;
    const int numPiMinus = 2;
    const int rhoPrimePDG = 30113;

    // Variables for the output tree
    float mass = -1.0f;
    float ptVal = -1.0f;
    float rapidity = -10.0f;
    float phiVal = -10.0f;

    int trackSigns[4] = {0, 0, 0, 0};
    float trackPts[4] = {-1.0f, -1.0f, -1.0f, -1.0f};
    float trackEtas[4] = {-10.0f, -10.0f, -10.0f, -10.0f};
    float trackPhis[4] = {-10.0f, -10.0f, -10.0f, -10.0f};
    int trackPdgs[4] = {0, 0, 0, 0};
    bool isPrimary[4] = {false, false, false, false};

    // Count all MC events
    mcRegistry.fill(HIST("MC/Events/hAllEvents"), 0);

    bool foundValidRhoPrime = false;
    bool passedCuts = false;

    // Search for rho prime particles in MC
    for (const auto& particle : mcParticles) {
      if (particle.pdgCode() != rhoPrimePDG) {
        continue;
      }

      // Get daughters of the rho prime candidate
      auto daughters = particle.daughters_as<UDMcParticles>();
      mcRegistry.fill(HIST("MC/Control/hNDaughters"), daughters.size());

      // Require exactly 4 daughters
      if (daughters.size() != numFourPionTracks) {
        continue;
      }

      // Reset arrays for new candidate
      for (int i = 0; i < 4; i++) {
        trackSigns[i] = 0;
        trackPts[i] = -1.0f;
        trackEtas[i] = -10.0f;
        trackPhis[i] = -10.0f;
        trackPdgs[i] = 0;
        isPrimary[i] = false;
      }

      // Reconstruct the 4-pion system from daughters
      PxPyPzMVector fourPionSystem;
      int nPiPlus = 0, nPiMinus = 0;
      int pionIndex = 0;

      for (const auto& daughter : daughters) {
        PxPyPzMVector pionVector(daughter.px(), daughter.py(), daughter.pz(),
                                 o2::constants::physics::MassPionCharged);

        if (daughter.pdgCode() == 211) { // pi plus
          nPiPlus++;

          // Fill arrays for the tree
          if (pionIndex < 4) {
            trackSigns[pionIndex] = 1;
            trackPts[pionIndex] = pt(daughter.px(), daughter.py());
            trackEtas[pionIndex] = eta(daughter.px(), daughter.py(), daughter.pz());
            trackPhis[pionIndex] = phi(daughter.px(), daughter.py());
            trackPdgs[pionIndex] = daughter.pdgCode();
            isPrimary[pionIndex] = daughter.isPhysicalPrimary();
            pionIndex++;
          }

          fourPionSystem += pionVector;

          // Fill individual pion histograms
          mcRegistry.fill(HIST("MC/Tracks/hPdgCode"), daughter.pdgCode());
          mcRegistry.fill(HIST("MC/Tracks/hPt"), trackPts[pionIndex - 1]);
          mcRegistry.fill(HIST("MC/Tracks/hEta"), trackEtas[pionIndex - 1]);
          mcRegistry.fill(HIST("MC/Tracks/hPhi"), trackPhis[pionIndex - 1]);

        } else if (daughter.pdgCode() == -211) { // pi minus
          nPiMinus++;

          // Fill arrays for the tree
          if (pionIndex < 4) {
            trackSigns[pionIndex] = -1;
            trackPts[pionIndex] = pt(daughter.px(), daughter.py());
            trackEtas[pionIndex] = eta(daughter.px(), daughter.py(), daughter.pz());
            trackPhis[pionIndex] = phi(daughter.px(), daughter.py());
            trackPdgs[pionIndex] = daughter.pdgCode();
            isPrimary[pionIndex] = daughter.isPhysicalPrimary();
            pionIndex++;
          }

          fourPionSystem += pionVector;
          // Fill individual pion histograms
          mcRegistry.fill(HIST("MC/Tracks/hPdgCode"), daughter.pdgCode());
          mcRegistry.fill(HIST("MC/Tracks/hPt"), trackPts[pionIndex - 1]);
          mcRegistry.fill(HIST("MC/Tracks/hEta"), trackEtas[pionIndex - 1]);
          mcRegistry.fill(HIST("MC/Tracks/hPhi"), trackPhis[pionIndex - 1]);
        }
      }

      // Fill control histograms for pion composition
      mcRegistry.fill(HIST("MC/Control/hNPiPlus"), nPiPlus);
      mcRegistry.fill(HIST("MC/Control/hNPiMinus"), nPiMinus);

      // Require exactly 2pi plus and 2pi minus for neutral system
      if (nPiPlus != numPiPlus || nPiMinus != numPiMinus) {
        continue;
      }

      // Calculate system properties
      mass = fourPionSystem.M();
      ptVal = fourPionSystem.Pt();
      rapidity = fourPionSystem.Rapidity();
      phiVal = fourPionSystem.Phi();

      // Fill control histograms
      mcRegistry.fill(HIST("MC/Control/hMAll"), mass);
      mcRegistry.fill(HIST("MC/Control/hPtAll"), ptVal);
      mcRegistry.fill(HIST("MC/Control/hYAll"), rapidity);
      mcRegistry.fill(HIST("MC/Control/hMvsPtAll"), mass, ptVal);

      // Apply kinematic cuts
      passedCuts = true;

      if (mass < systemMassMinCut || mass > systemMassMaxCut) {
        mcRegistry.fill(HIST("MC/Cuts/hMassRejected"), mass);
        passedCuts = false;
      }
      if (ptVal > systemPtCut) {
        mcRegistry.fill(HIST("MC/Cuts/hPtRejected"), ptVal);
        passedCuts = false;
      }
      if (std::abs(rapidity) > systemYCut) {
        mcRegistry.fill(HIST("MC/Cuts/hYRejected"), rapidity);
        passedCuts = false;
      }

      foundValidRhoPrime = true;

      // Fill final histograms only for events passing all cuts
      if (passedCuts) {
        mcRegistry.fill(HIST("MC/System/hM"), mass);
        mcRegistry.fill(HIST("MC/System/hPt"), ptVal);
        mcRegistry.fill(HIST("MC/System/hY"), rapidity);
        mcRegistry.fill(HIST("MC/System/hMvsPt"), mass, ptVal);
        mcRegistry.fill(HIST("MC/System/hMvsY"), mass, rapidity);
        mcRegistry.fill(HIST("MC/Events/hAccepted"), 0);
      }

      // Fill rho prime specific histograms
      mcRegistry.fill(HIST("MC/Events/hNPions"), 4);
      mcRegistry.fill(HIST("MC/RhoPrime/hFound"), 1);
      mcRegistry.fill(HIST("MC/RhoPrime/hDecayPions"), 4);
      mcRegistry.fill(HIST("MC/RhoPrime/hMass"), mass);
      mcRegistry.fill(HIST("MC/RhoPrime/hPt"), ptVal);

      break; // Process only the first rho prime found per event
    }

    // Fill output tree only for valid rho prime events passing all cuts
    if (foundValidRhoPrime && passedCuts) {
      mcFourPiTree(
        runIndex,
        mass, ptVal, rapidity, phiVal,
        0.0f, 0.0f, 0.0f, // Vertex position
        trackSigns, trackPts, trackEtas, trackPhis,
        trackPdgs,
        isPrimary[0], isPrimary[1], isPrimary[2], isPrimary[3]);
    }
  }

  PROCESS_SWITCH(upcRhoPrimeAnalysis, processRealData, "Process real data", true);
  PROCESS_SWITCH(upcRhoPrimeAnalysis, processMCgen, "Process MC generated data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<upcRhoPrimeAnalysis>(cfgc, TaskName{"upc-rho-prime-analysis"})};
}
