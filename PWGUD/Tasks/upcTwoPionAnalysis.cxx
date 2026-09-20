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
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(M, m, double);
DECLARE_SOA_COLUMN(Pt, pt, double);
DECLARE_SOA_COLUMN(Y, y, double);
DECLARE_SOA_COLUMN(Eta, eta, double);
DECLARE_SOA_COLUMN(Phi, phi, double);
DECLARE_SOA_COLUMN(PosX, posX, double);
DECLARE_SOA_COLUMN(PosY, posY, double);
DECLARE_SOA_COLUMN(PosZ, posZ, double);
DECLARE_SOA_COLUMN(NumContrib, numContrib, int32_t);
DECLARE_SOA_COLUMN(TotalFT0AmplitudeA, totalFT0AmplitudeA, float);
DECLARE_SOA_COLUMN(TotalFT0AmplitudeC, totalFT0AmplitudeC, float);
DECLARE_SOA_COLUMN(TotalFV0AmplitudeA, totalFV0AmplitudeA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNA, energyCommonZNA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNC, energyCommonZNC, float);
DECLARE_SOA_COLUMN(Sign, sign, std::vector<int>);
DECLARE_SOA_COLUMN(TrackPt, trackPt, std::vector<float>);
DECLARE_SOA_COLUMN(TrackEta, trackEta, std::vector<float>);
DECLARE_SOA_COLUMN(TrackPhi, trackPhi, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNSigmaEl, tpcNSigmaEl, std::vector<float>);
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, std::vector<float>);
DECLARE_SOA_COLUMN(IsReconstructedWithUPC, isReconstructedWithUPC, bool);
} // namespace twopi

// Define the output
DECLARE_SOA_TABLE(SYSTEMTREE, "AOD", "SystemTree",
                  twopi::RunNumber, twopi::M, twopi::Pt, twopi::Y, twopi::Eta, twopi::Phi,
                  twopi::PosX, twopi::PosY, twopi::PosZ, twopi::NumContrib,
                  twopi::TotalFT0AmplitudeA, twopi::TotalFT0AmplitudeC, twopi::TotalFV0AmplitudeA,
                  twopi::EnergyCommonZNA, twopi::EnergyCommonZNC,
                  twopi::Sign, twopi::TrackPt, twopi::TrackEta, twopi::TrackPhi,
                  twopi::TPCNSigmaEl, twopi::TPCNSigmaPi, twopi::IsReconstructedWithUPC)
} // namespace o2::aod

struct upcTwoPionAnalysis {
  Produces<aod::SYSTEMTREE> systemTree;

  // ==================== SYSTEM SELECTION ====================
  Configurable<double> systemYCut{"systemYCut", 0.9, "Rapidity cut for reco system"};
  Configurable<double> systemMassMinCut{"systemMassMinCut", 0.5, "Min M cut for reco system"};
  Configurable<double> systemMassMaxCut{"systemMassMaxCut", 1.0, "Max M cut for reco system"};
  Configurable<double> systemPtCut{"systemPtCut", 0.1, "Max pT cut for reco system"};

  // ==================== EVENT SELECTION ====================
  Configurable<float> vZCut{"vZCut", 10.0, "Cut on vertex Z position"};
  Configurable<int> numPVContrib{"numPVContrib", 2, "Number of PV contributors"};
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
  Configurable<double> etaCut{"etaCut", 0.9, "Track Pseudorapidity"};

  // ==================== GAP SELECTION ====================
  Configurable<bool> cutGapSide{"cutGapSide", true, "apply gap side cut?"};
  Configurable<bool> useTrueGap{"useTrueGap", false, "use true gap?"};
  Configurable<float> cutTrueGapSideFV0{"cutTrueGapSideFV0", 180000, "FV0A threshold for SG selector"};
  Configurable<float> cutTrueGapSideFT0A{"cutTrueGapSideFT0A", 150., "FT0A threshold for SG selector"};
  Configurable<float> cutTrueGapSideFT0C{"cutTrueGapSideFT0C", 50., "FT0C threshold for SG selector"};
  Configurable<float> cutTrueGapSideZDC{"cutTrueGapSideZDC", 10000., "ZDC threshold for SG selector"};

  // ==================== RCT FLAGS ====================
  Configurable<bool> useRctFlag{"useRctFlag", true, "use RCT flags for event selection?"};
  Configurable<int> cutRctFlag{"cutRctFlag", 1, "0=off, 1=CBT, 2=CBT+ZDC"};

  // ==================== TRACK SELECTION ====================
  Configurable<bool> useOnlyPVtracks{"useOnlyPVtracks", true, "Use only PV tracks"};
  Configurable<float> tpcChi2NClsCut{"tpcChi2NClsCut", 4.0, "TPC chi2/N clusters cut"};
  Configurable<float> itsChi2NClsCut{"itsChi2NClsCut", 36.0, "ITS chi2/N clusters cut"};
  Configurable<float> nSigmaTPCcut{"nSigmaTPCcut", 3.0, "TPC nSigma cut"};
  Configurable<float> dcaXYcut{"dcaXYcut", 0.0105, "dcaXY cut"};
  Configurable<float> dcaZcut{"dcaZcut", 2.0, "dcaZ cut"};
  Configurable<int> minTPCFindableClusters{"minTPCFindableClusters", 70, "Minimum number of findable TPC clusters"};
  Configurable<int> tracksMinTpcNClsCrossedRowsCut{"tracksMinTpcNClsCrossedRowsCut", 130, "min TPC crossed rows"};
  Configurable<int> tracksMinTpcNClsCut{"tracksMinTpcNClsCut", 120, "min TPC clusters"};
  Configurable<int> tracksMinItsNClsCut{"tracksMinItsNClsCut", 4, "min ITS clusters"};
  Configurable<float> tracksMinTpcChi2NClCut{"tracksMinTpcChi2NClCut", 1.0, "min TPC chi2/Ncls"};
  Configurable<float> tracksMaxTpcChi2NClCut{"tracksMaxTpcChi2NClCut", 3.0, "max TPC chi2/Ncls"};
  Configurable<float> tracksDcaMaxCut{"tracksDcaMaxCut", 1.0, "max DCA cut on tracks"};
  Configurable<float> tracksMinPtCut{"tracksMinPtCut", 0.1, "min pT cut"};
  Configurable<float> tracksMinTpcNClsCrossedOverFindableCut{"tracksMinTpcNClsCrossedOverFindableCut", 1.0, "min crossed/findable"};

  // ==================== PID CONFIGURATION ====================
  Configurable<float> tracksTpcNSigmaPiCut{"tracksTpcNSigmaPiCut", 3.0, "TPC nSigma pion cut"};
  Configurable<bool> rejectLowerProbPairs{"rejectLowerProbPairs", true, "reject track pairs with smaller PID radii"};
  Configurable<bool> requireTof{"requireTof", false, "require TOF signal?"};

  // ==================== HISTOGRAMS ====================
  HistogramRegistry registry{
    "registry",
    {
      {"Events/Flow", "Event flow;Cut;Counts", {HistType::kTH1F, {{12, 0, 12}}}},
      {"Events/VertexZ", "Vertex Z;z (cm);Counts", {HistType::kTH1F, {{200, -20, 20}}}},
      {"Events/NumContrib", "Number of contributors;N_{contrib};Counts", {HistType::kTH1F, {{100, 0, 100}}}},
      {"Events/FV0Amplitude", "FV0 amplitude;Amplitude;Counts", {HistType::kTH1F, {{200, 0, 200}}}},
      {"Events/FT0AmplitudeA", "FT0A amplitude;Amplitude;Counts", {HistType::kTH1F, {{200, 0, 200}}}},
      {"Events/FT0AmplitudeC", "FT0C amplitude;Amplitude;Counts", {HistType::kTH1F, {{200, 0, 200}}}},
      {"Events/ZDCEnergy", "ZDC energy;Energy (TeV);Counts", {HistType::kTH1F, {{200, 0, 2}}}},

      {"Tracks/Pt", "Track p_{T};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{200, 0, 1.5}}}},
      {"Tracks/Eta", "Track #eta;#eta;Counts", {HistType::kTH1F, {{200, -2, 2}}}},
      {"Tracks/TPCNSigmaPi", "TPC n#sigma for #pi;n#sigma;Counts", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Tracks/TPCChi2NCl", "TPC #chi^{2}/N_{cls};#chi^{2}/N_{cls};Counts", {HistType::kTH1F, {{200, 0, 20}}}},
      {"Tracks/ITSChi2NCl", "ITS #chi^{2}/N_{cls};#chi^{2}/N_{cls};Counts", {HistType::kTH1F, {{200, 0, 50}}}},
      {"Tracks/RejectionReasons", "Track rejection reasons;Reason;Counts", {HistType::kTH1F, {{15, 0, 15}}}},
      {"Tracks/DCASpectrum", "Track DCA spectrum;DCA (cm);Counts", {HistType::kTH1F, {{100, 0, 5}}}},
      {"Tracks/ChargeDistribution", "Track charge distribution;Charge;Counts", {HistType::kTH1F, {{3, -1.5, 1.5}}}},
      {"Tracks/TPCClusters", "TPC clusters findable;N_{clusters};Counts", {HistType::kTH1F, {{100, 0, 200}}}},

      {"System/hM", ";m (GeV/#it{c}^{2});counts", {HistType::kTH1F, {{1000, 0.0, 2.5}}}},
      {"System/hPt", ";p_{T} (GeV/#it{c});counts", {HistType::kTH1F, {{1000, 0.0, 1.0}}}},
      {"System/hEta", ";#eta;counts", {HistType::kTH1F, {{180, -1.0, 1.0}}}},
      {"System/hPhi", ";#phi;counts", {HistType::kTH1F, {{180, 0.0, 6.28}}}},
      {"System/hY", ";y;counts", {HistType::kTH1F, {{180, -1.0, 1.0}}}},

      {"Cuts/MBefore", "Mass before cuts;m (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0, 2.5}}}},
      {"Cuts/MAfter", "Mass after cuts;m (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0, 2.5}}}},
      {"Cuts/PtBefore", "p_{T} before cuts;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0, 1.0}}}},
      {"Cuts/PtAfter", "p_{T} after cuts;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0, 1.0}}}}
    }
  };

  void init(InitContext&)
  {
    auto hFlow = registry.get<TH1>(HIST("Events/Flow"));
    hFlow->GetXaxis()->SetBinLabel(1, "All events");
    hFlow->GetXaxis()->SetBinLabel(2, "ITS-TPC cut");
    hFlow->GetXaxis()->SetBinLabel(3, "SBP cut");
    hFlow->GetXaxis()->SetBinLabel(4, "ITS ROFb cut");
    hFlow->GetXaxis()->SetBinLabel(5, "TFB cut");
    hFlow->GetXaxis()->SetBinLabel(6, "Gap Side cut");
    hFlow->GetXaxis()->SetBinLabel(7, "ZDC cut");
    hFlow->GetXaxis()->SetBinLabel(8, "FV0/FT0 cut");
    hFlow->GetXaxis()->SetBinLabel(9, "PV contrib cut");
    hFlow->GetXaxis()->SetBinLabel(10, "Z vtx cut");
    hFlow->GetXaxis()->SetBinLabel(11, "2 tracks found");
    hFlow->GetXaxis()->SetBinLabel(12, "System cuts");

    auto hReject = registry.get<TH1>(HIST("Tracks/RejectionReasons"));
    hReject->GetXaxis()->SetBinLabel(1, "All Tracks");
    hReject->GetXaxis()->SetBinLabel(2, "PV Contributor");
    hReject->GetXaxis()->SetBinLabel(3, "Has ITS+TPC");
    hReject->GetXaxis()->SetBinLabel(4, "pT cut");
    hReject->GetXaxis()->SetBinLabel(5, "TPC chi2/cluster");
    hReject->GetXaxis()->SetBinLabel(6, "ITS chi2/cluster");
    hReject->GetXaxis()->SetBinLabel(7, "TPC crossed rows");
    hReject->GetXaxis()->SetBinLabel(8, "TPC clusters");
    hReject->GetXaxis()->SetBinLabel(9, "TPC crossed/findable");
    hReject->GetXaxis()->SetBinLabel(10, "TPC nSigmaPi");
    hReject->GetXaxis()->SetBinLabel(11, "Eta acceptance");
    hReject->GetXaxis()->SetBinLabel(12, "DCAz cut");
    hReject->GetXaxis()->SetBinLabel(13, "DCAxy cut");
    hReject->GetXaxis()->SetBinLabel(14, "Accepted Tracks");
  }

  static double eta(double px, double py, double pz)
  {
    double p = std::sqrt(px * px + py * py + pz * pz);
    return 0.5 * std::log((p + pz) / (p - pz + 1e-12));
  }

  static double phi(double px, double py)
  {
    return std::atan2(py, px);
  }

  void process(UDCollisions::iterator const& collision, UDtracks const& tracks)
  {
    registry.fill(HIST("Events/Flow"), 0);

    registry.fill(HIST("Events/VertexZ"), collision.posZ());
    registry.fill(HIST("Events/NumContrib"), collision.numContrib());
    registry.fill(HIST("Events/FV0Amplitude"), collision.totalFV0AmplitudeA());
    registry.fill(HIST("Events/FT0AmplitudeA"), collision.totalFT0AmplitudeA());
    registry.fill(HIST("Events/FT0AmplitudeC"), collision.totalFT0AmplitudeC());
    registry.fill(HIST("Events/ZDCEnergy"), collision.energyCommonZNA());
    registry.fill(HIST("Events/ZDCEnergy"), collision.energyCommonZNC());

    // ==================== EVENT SELECTION ====================
    if (!collision.vtxITSTPC())
      return;
    registry.fill(HIST("Events/Flow"), 1);

    if (!collision.sbp())
      return;
    registry.fill(HIST("Events/Flow"), 2);

    if (!collision.itsROFb())
      return;
    registry.fill(HIST("Events/Flow"), 3);

    if (!collision.tfb())
      return;
    registry.fill(HIST("Events/Flow"), 4);

    if (collision.gapSide() != gapSide)
      return;
    registry.fill(HIST("Events/Flow"), 5);

    // ZDC cut
    if (collision.energyCommonZNA() > zdcCut || collision.energyCommonZNC() > zdcCut)
      return;
    registry.fill(HIST("Events/Flow"), 6);

    // FIT cuts
    if (collision.totalFV0AmplitudeA() > fv0Cut ||
        collision.totalFT0AmplitudeA() > ft0aCut ||
        collision.totalFT0AmplitudeC() > ft0cCut)
      return;
    registry.fill(HIST("Events/Flow"), 7);

    if (collision.numContrib() != numPVContrib)
      return;
    registry.fill(HIST("Events/Flow"), 8);

    if (std::abs(collision.posZ()) > vZCut)
      return;
    registry.fill(HIST("Events/Flow"), 9);

    // ==================== TRACK SELECTION ====================
    std::vector<decltype(tracks.begin())> posPions;
    std::vector<decltype(tracks.begin())> negPions;
    posPions.reserve(2);
    negPions.reserve(2);

    for (const auto& track : tracks) {
      registry.fill(HIST("Tracks/RejectionReasons"), 0);

      // PV contributor
      if (useOnlyPVtracks && !track.isPVContributor()) {
        registry.fill(HIST("Tracks/RejectionReasons"), 1);
        continue;
      }

      // Has ITS+TPC
      if (!track.hasITS() || !track.hasTPC()) {
        registry.fill(HIST("Tracks/RejectionReasons"), 2);
        continue;
      }

      // Fill diagnostic histograms
      registry.fill(HIST("Tracks/Pt"), track.pt());
      registry.fill(HIST("Tracks/Eta"), eta(track.px(), track.py(), track.pz()));
      registry.fill(HIST("Tracks/TPCNSigmaPi"), track.tpcNSigmaPi());
      registry.fill(HIST("Tracks/TPCChi2NCl"), track.tpcChi2NCl());
      registry.fill(HIST("Tracks/ITSChi2NCl"), track.itsChi2NCl());
      registry.fill(HIST("Tracks/DCASpectrum"), std::hypot(track.dcaXY(), track.dcaZ()));
      registry.fill(HIST("Tracks/ChargeDistribution"), track.sign());
      registry.fill(HIST("Tracks/TPCClusters"), track.tpcNClsFindable());

      // pT cut
      if (track.pt() <= tracksMinPtCut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 3);
        continue;
      }

      // TPC chi2/Ncls
      if (track.tpcChi2NCl() < tracksMinTpcChi2NClCut ||
          track.tpcChi2NCl() > tracksMaxTpcChi2NClCut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 4);
        continue;
      }

      // ITS chi2/Ncls
      if (track.itsChi2NCl() > itsChi2NClsCut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 5);
        continue;
      }

      // TPC crossed rows
      if (track.tpcNClsCrossedRows() < tracksMinTpcNClsCrossedRowsCut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 6);
        continue;
      }

      // TPC clusters
      if (track.tpcNClsFindable() < tracksMinTpcNClsCut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 7);
        continue;
      }

      // TPC crossed/findable
      float crossedOverFindable = static_cast<float>(track.tpcNClsCrossedRows()) /
                                   std::max<int>(1, track.tpcNClsFindable());
      if (crossedOverFindable < tracksMinTpcNClsCrossedOverFindableCut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 8);
        continue;
      }

      // TPC nSigma for pions
      if (std::abs(track.tpcNSigmaPi()) > nSigmaTPCcut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 9);
        continue;
      }

      // Eta acceptance
      float trackEta = eta(track.px(), track.py(), track.pz());
      if (std::abs(trackEta) > etaCut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 10);
        continue;
      }

      // DCAz
      if (std::abs(track.dcaZ()) > dcaZcut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 11);
        continue;
      }

      // DCAxy - paramétrico
      float maxDCAxy = 0.0105 + 0.035 / std::pow(track.pt(), 1.1);
      if (std::abs(track.dcaXY()) > maxDCAxy) {
        registry.fill(HIST("Tracks/RejectionReasons"), 12);
        continue;
      }

      registry.fill(HIST("Tracks/RejectionReasons"), 13);

      // Store tracks
      if (track.sign() > 0 && posPions.size() < 2) {
        posPions.push_back(track);
      } else if (track.sign() < 0 && negPions.size() < 2) {
        negPions.push_back(track);
      }
    }

    if (posPions.size() < 1 || negPions.size() < 1)
      return;
    registry.fill(HIST("Events/Flow"), 10);

    // ==================== SYSTEM RECONSTRUCTION ====================
    ROOT::Math::PxPyPzMVector twoPionSystem;

    for (const auto& track : {posPions[0], negPions[0]}) {
      ROOT::Math::PxPyPzMVector pionVec(
        track.px(), track.py(), track.pz(),
        o2::constants::physics::MassPionCharged);
      twoPionSystem += pionVec;
    }

    registry.fill(HIST("Cuts/MBefore"), twoPionSystem.M());
    registry.fill(HIST("Cuts/PtBefore"), twoPionSystem.Pt());

    // ==================== SYSTEM CUTS ====================
    if (twoPionSystem.M() < systemMassMinCut || twoPionSystem.M() > systemMassMaxCut)
      return;
    if (std::abs(twoPionSystem.Rapidity()) > systemYCut)
      return;
    if (twoPionSystem.Pt() > systemPtCut)
      return;
    registry.fill(HIST("Events/Flow"), 11);

    registry.fill(HIST("Cuts/MAfter"), twoPionSystem.M());
    registry.fill(HIST("Cuts/PtAfter"), twoPionSystem.Pt());
    registry.fill(HIST("System/hM"), twoPionSystem.M());
    registry.fill(HIST("System/hPt"), twoPionSystem.Pt());
    registry.fill(HIST("System/hEta"), twoPionSystem.Eta());
    registry.fill(HIST("System/hPhi"), twoPionSystem.Phi() + o2::constants::math::PI);
    registry.fill(HIST("System/hY"), twoPionSystem.Rapidity());

    // ==================== FILL OUTPUT ====================
    std::vector<float> trackPts, trackEtas, trackPhis;
    std::vector<int> trackSigns;
    std::vector<float> tpcNSigmasEl, tpcNSigmasPi;

    // Store info for both selected tracks
    for (const auto& track : {posPions[0], negPions[0]}) {
      trackPts.push_back(track.pt());
      trackEtas.push_back(eta(track.px(), track.py(), track.pz()));
      trackPhis.push_back(phi(track.px(), track.py()));
      trackSigns.push_back(track.sign());
      tpcNSigmasEl.push_back(track.tpcNSigmaEl());
      tpcNSigmasPi.push_back(track.tpcNSigmaPi());
    }

    bool isReconstructedWithUPC = (collision.flags() == 1);

    systemTree(
      collision.runNumber(),
      twoPionSystem.M(),
      twoPionSystem.Pt(),
      twoPionSystem.Rapidity(),
      twoPionSystem.Eta(),
      twoPionSystem.Phi(),
      collision.posX(),
      collision.posY(),
      collision.posZ(),
      collision.numContrib(),
      collision.totalFT0AmplitudeA(),
      collision.totalFT0AmplitudeC(),
      collision.totalFV0AmplitudeA(),
      collision.energyCommonZNA(),
      collision.energyCommonZNC(),
      trackSigns,
      trackPts,
      trackEtas,
      trackPhis,
      tpcNSigmasEl,
      tpcNSigmasPi,
      isReconstructedWithUPC
    );
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<upcTwoPionAnalysis>(cfgc, TaskName{"upc-two-pion-analysis"})};
}