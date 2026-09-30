///
/// \file   upcTwoPionAnalysis.cxx
/// \brief  Analisis basico de fotoproduccion coherente de dipiones (rho0) en UPC
///         usando tablas UD derivadas del SG-producer.
///         En el arbol solo se guardan pares pi+ pi- que son los DOS UNICOS
///         contribuyentes del mismo vertice primario.
///

#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"
#include "PWGUD/DataModel/UDTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h>
#include <TH1.h>
#include <TH2.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::framework;

using UDtracks = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
using UDCollisions = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDZdcsReduced>;

namespace o2::aod
{
namespace twopi
{
// ---- evento ----
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(RecoFlag, recoFlag, int); // 0 = std, 1 = UPC
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(NumContrib, numContrib, int);
DECLARE_SOA_COLUMN(TotalFT0AmplitudeA, totalFT0AmplitudeA, float);
DECLARE_SOA_COLUMN(TotalFT0AmplitudeC, totalFT0AmplitudeC, float);
DECLARE_SOA_COLUMN(TotalFV0AmplitudeA, totalFV0AmplitudeA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNA, energyCommonZNA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNC, energyCommonZNC, float);
DECLARE_SOA_COLUMN(TimeZNA, timeZNA, float);
DECLARE_SOA_COLUMN(TimeZNC, timeZNC, float);
DECLARE_SOA_COLUMN(NeutronClass, neutronClass, int); // 0 = 0n0n, 1 = Xn0n, 2 = 0nXn, 3 = XnXn
// ---- sistema pi+ pi- ----
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Phi, phi, float);              // [0, 2pi), igual que en el histograma
DECLARE_SOA_COLUMN(IsSelected, isSelected, bool); // pasa los cortes finales en m, pT, y
// ---- pi+ ----
DECLARE_SOA_COLUMN(PiPlusPt, piPlusPt, float);
DECLARE_SOA_COLUMN(PiPlusEta, piPlusEta, float);
DECLARE_SOA_COLUMN(PiPlusPhi, piPlusPhi, float);
DECLARE_SOA_COLUMN(PiPlusTpcNSigmaPi, piPlusTpcNSigmaPi, float);
DECLARE_SOA_COLUMN(PiPlusTpcNSigmaEl, piPlusTpcNSigmaEl, float);
// ---- pi- ----
DECLARE_SOA_COLUMN(PiMinusPt, piMinusPt, float);
DECLARE_SOA_COLUMN(PiMinusEta, piMinusEta, float);
DECLARE_SOA_COLUMN(PiMinusPhi, piMinusPhi, float);
DECLARE_SOA_COLUMN(PiMinusTpcNSigmaPi, piMinusTpcNSigmaPi, float);
DECLARE_SOA_COLUMN(PiMinusTpcNSigmaEl, piMinusTpcNSigmaEl, float);
} // namespace twopi

DECLARE_SOA_TABLE(SYSTEMTREE, "AOD", "SystemTree",
                  twopi::RunNumber, twopi::RecoFlag,
                  twopi::PosX, twopi::PosY, twopi::PosZ, twopi::NumContrib,
                  twopi::TotalFT0AmplitudeA, twopi::TotalFT0AmplitudeC, twopi::TotalFV0AmplitudeA,
                  twopi::EnergyCommonZNA, twopi::EnergyCommonZNC, twopi::TimeZNA, twopi::TimeZNC, twopi::NeutronClass,
                  twopi::M, twopi::Pt, twopi::Y, twopi::Phi, twopi::IsSelected,
                  twopi::PiPlusPt, twopi::PiPlusEta, twopi::PiPlusPhi, twopi::PiPlusTpcNSigmaPi, twopi::PiPlusTpcNSigmaEl,
                  twopi::PiMinusPt, twopi::PiMinusEta, twopi::PiMinusPhi, twopi::PiMinusTpcNSigmaPi, twopi::PiMinusTpcNSigmaEl);
} // namespace o2::aod

struct UpcTwoPionAnalysis {
  Produces<aod::SYSTEMTREE> systemTree;
  SGSelector sgSelector;

  // ==================== EVENT SELECTION ====================
  Configurable<bool> cutGapSide{"cutGapSide", true, "apply gap side cut?"};
  Configurable<int> gapSide{"gapSide", 2, "required gap side (0 = A, 1 = C, 2 = double gap)"};
  Configurable<bool> useTrueGap{"useTrueGap", false, "recompute gap with SGSelector::trueGap?"};
  Configurable<float> cutTrueGapSideFV0{"cutTrueGapSideFV0", 180000, "FV0A threshold for trueGap"};
  Configurable<float> cutTrueGapSideFT0A{"cutTrueGapSideFT0A", 150., "FT0A threshold for trueGap"};
  Configurable<float> cutTrueGapSideFT0C{"cutTrueGapSideFT0C", 50., "FT0C threshold for trueGap"};
  Configurable<float> cutTrueGapSideZDC{"cutTrueGapSideZDC", 10000., "ZDC threshold for trueGap"};

  Configurable<bool> vtxITSTPCcut{"vtxITSTPCcut", true, "require ITS-TPC vertex"};
  Configurable<bool> sbpCut{"sbpCut", true, "reject same-bunch pile-up"};
  Configurable<bool> itsROFbCut{"itsROFbCut", true, "reject ITS ROF border"};
  Configurable<bool> tfbCut{"tfbCut", true, "reject TF border"};

  Configurable<bool> useRctFlag{"useRctFlag", true, "use RCT flags?"};
  Configurable<int> cutRctFlag{"cutRctFlag", 1, "0 = off, 1 = CBT, 2 = CBT+ZDC"};

  Configurable<float> vZCut{"vZCut", 10.0, "max |z_vtx| (cm)"};
  Configurable<int> numPVContrib{"numPVContrib", 2, "required number of PV contributors"};

  Configurable<bool> cutFit{"cutFit", true, "apply FIT amplitude veto?"};
  Configurable<float> fv0Cut{"fv0Cut", 50.0, "max FV0A amplitude"};
  Configurable<float> ft0aCut{"ft0aCut", 50.0, "max FT0A amplitude"};
  Configurable<float> ft0cCut{"ft0cCut", 50.0, "max FT0C amplitude"};

  Configurable<float> znTimeCut{"znTimeCut", 2.0, "|ZN time| (ns) below which a ZN is considered hit"};
  Configurable<int> neutronClassCut{"neutronClassCut", 0, "-1 = all, 0 = 0n0n, 1 = Xn0n, 2 = 0nXn, 3 = XnXn"};

  // ==================== TRACK SELECTION ====================
  // Nota: isPVContributor() es OBLIGATORIO (no configurable) -> ambos piones del mismo vertice
  Configurable<float> tracksMinPtCut{"tracksMinPtCut", 0.1, "min track pT (GeV/c)"};
  Configurable<float> etaCut{"etaCut", 0.9, "max track |eta|"};
  Configurable<int> tracksMinItsNClsCut{"tracksMinItsNClsCut", 4, "min ITS clusters"};
  Configurable<float> itsChi2NClsCut{"itsChi2NClsCut", 3.0, "max ITS chi2/Ncls"};
  Configurable<int> tracksMinTpcNClsCut{"tracksMinTpcNClsCut", 120, "min FOUND TPC clusters"};
  Configurable<int> tracksMinTpcNClsCrossedRowsCut{"tracksMinTpcNClsCrossedRowsCut", 130, "min TPC crossed rows"};
  Configurable<float> tracksMinTpcChi2NClCut{"tracksMinTpcChi2NClCut", 1.0, "min TPC chi2/Ncls"};
  Configurable<float> tracksMaxTpcChi2NClCut{"tracksMaxTpcChi2NClCut", 3.0, "max TPC chi2/Ncls"};
  Configurable<float> tracksMinTpcNClsCrossedOverFindableCut{"tracksMinTpcNClsCrossedOverFindableCut", 1.0, "min crossed rows / findable"};
  Configurable<float> dcaZcut{"dcaZcut", 1.0, "max |DCAz| (cm)"};
  Configurable<bool> requireTof{"requireTof", false, "require TOF signal?"};

  // ==================== PID ====================
  Configurable<float> tracksTpcNSigmaPiCut{"tracksTpcNSigmaPiCut", 3.0, "max pair PID radius for pions"};
  Configurable<bool> rejectLowerProbPairs{"rejectLowerProbPairs", true, "reject pairs closer to e/K/p than to pi"};

  // ==================== SYSTEM SELECTION ====================
  Configurable<float> systemMassMinCut{"systemMassMinCut", 0.5, "min m (GeV/c2)"};
  Configurable<float> systemMassMaxCut{"systemMassMaxCut", 1.0, "max m (GeV/c2)"};
  Configurable<float> systemPtCut{"systemPtCut", 0.1, "max pT (GeV/c): 0.1 PbPb, 0.3 OO/NeNe"};
  Configurable<float> systemYCut{"systemYCut", 0.9, "max |y|"};
  Configurable<bool> saveOnlySelected{"saveOnlySelected", false, "tree: true = solo candidatos que pasan m/pT/y; false = |y| y pT < treeMaxPtCut"};
  Configurable<float> treeMaxPtCut{"treeMaxPtCut", 1.0, "max pT stored in tree when saveOnlySelected = false"};

  // ==================== AXES ====================
  ConfigurableAxis mAxis{"mAxis", {250, 0.0, 2.5}, "#it{m}_{#pi#pi} (GeV/#it{c}^{2})"};
  ConfigurableAxis ptAxis{"ptAxis", {400, 0.0, 2.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis pt2Axis{"pt2Axis", {1000, 0.0, 1.0}, "#it{p}_{T}^{2} (GeV^{2}/#it{c}^{2})"};
  ConfigurableAxis yAxis{"yAxis", {120, -1.2, 1.2}, "#it{y}"};
  ConfigurableAxis etaAxis{"etaAxis", {120, -1.2, 1.2}, "#it{#eta}"};
  ConfigurableAxis phiAxis{"phiAxis", {180, 0.0, o2::constants::math::TwoPI}, "#it{#phi} (rad)"};
  ConfigurableAxis nSigmaAxis{"nSigmaAxis", {200, -10.0, 10.0}, "TPC #it{n#sigma}"};
  ConfigurableAxis znEnergyAxis{"znEnergyAxis", {250, -5.0, 20.0}, "ZN common energy (TeV)"};
  ConfigurableAxis znTimeAxis{"znTimeAxis", {200, -10.0, 10.0}, "ZN time (ns)"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr std::string_view Stage[2] = {"before/", "after/"};
  static constexpr std::string_view Selection[2] = {"all/", "selected/"};
  static constexpr std::string_view Charge[2] = {"unlikeSign/", "likeSign/"};

  void init(InitContext&)
  {
    // ---------------- EVENTOS ----------------
    const std::vector<std::string> evLabels = {"all", "gap side", "ITS-TPC vtx", "no SBP", "no ITS ROF border", "no TF border",
                                               "RCT", "|z_{vtx}|", "N_{PV contrib}", "FIT veto", "neutron class",
                                               "2 selected tracks", "#pi PID", "unlike-sign", "system cuts"};
    const int nEv = static_cast<int>(evLabels.size());
    registry.add("Events/hSelectionCounter", ";;events passing", kTH1D, {{nEv, -0.5, nEv - 0.5}});
    for (int i = 0; i < nEv; i++)
      registry.get<TH1>(HIST("Events/hSelectionCounter"))->GetXaxis()->SetBinLabel(i + 1, evLabels[i].c_str());

    registry.add("Events/before/hPosZ", ";#it{z}_{vtx} (cm);counts", kTH1D, {{200, -20.0, 20.0}});
    registry.add("Events/before/hPosXY", ";#it{x}_{vtx} (cm);#it{y}_{vtx} (cm);counts", kTH2D, {{200, -0.2, 0.2}, {200, -0.2, 0.2}});
    registry.add("Events/before/hNumContrib", ";#it{N}_{PV contrib};counts", kTH1D, {{51, -0.5, 50.5}});
    registry.add("Events/before/hFT0A", ";FT0A amplitude;counts", kTH1D, {{200, 0.0, 200.0}});
    registry.add("Events/before/hFT0C", ";FT0C amplitude;counts", kTH1D, {{200, 0.0, 200.0}});
    registry.add("Events/before/hFV0A", ";FV0A amplitude;counts", kTH1D, {{300, 0.0, 300.0}});
    registry.add("Events/before/hZNEnergy", ";ZNA common energy (TeV);ZNC common energy (TeV);counts", kTH2D, {znEnergyAxis, znEnergyAxis});
    registry.add("Events/before/hZNTime", ";ZNA time (ns);ZNC time (ns);counts", kTH2D, {znTimeAxis, znTimeAxis});
    registry.addClone("Events/before/", "Events/after/");

    // sanity check de "mismo vertice": deberia ser diagonal
    registry.add("Events/hNumContribVsPVTracks", ";PV-contributor tracks in UD table;collision.numContrib();counts", kTH2D, {{21, -0.5, 20.5}, {21, -0.5, 20.5}});
    registry.add("Events/hNeutronClass", ";;counts", kTH1D, {{4, -0.5, 3.5}});
    const char* ncLabels[4] = {"0n0n", "Xn0n", "0nXn", "XnXn"};
    for (int i = 0; i < 4; i++)
      registry.get<TH1>(HIST("Events/hNeutronClass"))->GetXaxis()->SetBinLabel(i + 1, ncLabels[i]);

    // ---------------- TRAZAS ----------------
    const std::vector<std::string> trkLabels = {"all", "PV contributor", "ITS hit", "ITS #it{N}_{cls}", "ITS IB hit", "ITS #chi^{2}/#it{N}_{cls}",
                                                "TPC hit", "TPC found #it{N}_{cls}", "TPC #chi^{2}/#it{N}_{cls}", "TPC crossed rows",
                                                "crossed/findable", "TOF", "#it{p}_{T}", "DCA", "#it{#eta}"};
    const int nTrk = static_cast<int>(trkLabels.size());
    registry.add("Tracks/hSelectionCounter", ";;tracks passing", kTH1D, {{nTrk, -0.5, nTrk - 0.5}});
    for (int i = 0; i < nTrk; i++)
      registry.get<TH1>(HIST("Tracks/hSelectionCounter"))->GetXaxis()->SetBinLabel(i + 1, trkLabels[i].c_str());

    registry.add("Tracks/before/hPt", ";#it{p}_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("Tracks/before/hEta", ";#it{#eta};counts", kTH1D, {etaAxis});
    registry.add("Tracks/before/hPhi", ";#it{#phi} (rad);counts", kTH1D, {phiAxis});
    registry.add("Tracks/before/hDcaXYVsPt", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);counts", kTH2D, {{100, 0.0, 2.0}, {200, -0.5, 0.5}});
    registry.add("Tracks/before/hDcaZ", ";DCA_{z} (cm);counts", kTH1D, {{200, -2.0, 2.0}});
    registry.add("Tracks/before/hItsNCls", ";ITS #it{N}_{cls};counts", kTH1D, {{8, -0.5, 7.5}});
    registry.add("Tracks/before/hItsChi2NCl", ";ITS #chi^{2}/#it{N}_{cls};counts", kTH1D, {{200, 0.0, 40.0}});
    registry.add("Tracks/before/hTpcNClsFound", ";TPC found #it{N}_{cls};counts", kTH1D, {{161, -0.5, 160.5}});
    registry.add("Tracks/before/hTpcCrossedRows", ";TPC crossed rows;counts", kTH1D, {{161, -0.5, 160.5}});
    registry.add("Tracks/before/hTpcChi2NCl", ";TPC #chi^{2}/#it{N}_{cls};counts", kTH1D, {{100, 0.0, 10.0}});
    registry.add("Tracks/before/hTpcSignalVsP", ";#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} (a.u.);counts", kTH2D, {{200, 0.0, 2.0}, {250, 0.0, 250.0}});
    registry.add("Tracks/before/hTpcNSigmaPi", ";TPC #it{n#sigma}_{#pi};counts", kTH1D, {nSigmaAxis});
    registry.addClone("Tracks/before/", "Tracks/after/");

    registry.add("Tracks/hNSelectedTracks", ";selected tracks per event;counts", kTH1D, {{11, -0.5, 10.5}});
    registry.add("Tracks/hTpcNSigmaPiVsEl", ";TPC #it{n#sigma}_{#pi};TPC #it{n#sigma}_{e};counts", kTH2D, {nSigmaAxis, nSigmaAxis});
    registry.add("Tracks/hPidRadiusPi", ";#it{n#sigma}_{#pi} pair radius;counts", kTH1D, {{200, 0.0, 20.0}});
    registry.add("Tracks/hPidRadiusEl", ";#it{n#sigma}_{e} pair radius;counts", kTH1D, {{200, 0.0, 20.0}});

    // ---------------- SISTEMA ----------------
    registry.add("System/all/unlikeSign/hM", ";#it{m}_{#pi#pi} (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("System/all/unlikeSign/hPt", ";#it{p}_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("System/all/unlikeSign/hPt2", ";#it{p}_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("System/all/unlikeSign/hY", ";#it{y};counts", kTH1D, {yAxis});
    registry.add("System/all/unlikeSign/hPhi", ";#it{#phi} (rad);counts", kTH1D, {phiAxis});
    registry.add("System/all/unlikeSign/hPtVsM", ";#it{m}_{#pi#pi} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.addClone("System/all/unlikeSign/", "System/all/likeSign/");
    registry.addClone("System/all/", "System/selected/");

    // N-1: cada variable con todos los cortes de sistema EXCEPTO el suyo (solo unlike-sign)
    registry.add("System/nMinus1/hM", ";#it{m}_{#pi#pi} (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("System/nMinus1/hPt", ";#it{p}_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("System/nMinus1/hPt2", ";#it{p}_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
  }

  // ---------- helpers ----------
  static double eta(double px, double py, double pz)
  {
    double p = std::sqrt(px * px + py * py + pz * pz);
    return 0.5 * std::log((p + pz) / (p - pz + 1e-12));
  }

  static double phi(double px, double py)
  {
    double ph = std::atan2(py, px);
    return ph < 0 ? ph + o2::constants::math::TwoPI : ph; // [0, 2pi)
  }

  int getNeutronClass(float timeZNA, float timeZNC) const
  {
    const bool hitA = std::abs(timeZNA) <= znTimeCut.value;
    const bool hitC = std::abs(timeZNC) <= znTimeCut.value;
    if (!hitA && !hitC)
      return 0; // 0n0n
    if (hitA && !hitC)
      return 1; // Xn0n
    if (!hitA && hitC)
      return 2; // 0nXn
    return 3;   // XnXn
  }

  template <typename C>
  bool isGoodRctFlag(const C& collision)
  {
    switch (cutRctFlag) {
      case 1:
        return sgSelector.isCBTOk(collision);
      case 2:
        return sgSelector.isCBTZdcOk(collision);
      default:
        return true;
    }
  }

  template <int stage, typename C>
  void fillEventQA(const C& collision, float znaE, float zncE, float znaT, float zncT)
  {
    registry.fill(HIST("Events/") + HIST(Stage[stage]) + HIST("hPosZ"), collision.posZ());
    registry.fill(HIST("Events/") + HIST(Stage[stage]) + HIST("hPosXY"), collision.posX(), collision.posY());
    registry.fill(HIST("Events/") + HIST(Stage[stage]) + HIST("hNumContrib"), collision.numContrib());
    registry.fill(HIST("Events/") + HIST(Stage[stage]) + HIST("hFT0A"), collision.totalFT0AmplitudeA());
    registry.fill(HIST("Events/") + HIST(Stage[stage]) + HIST("hFT0C"), collision.totalFT0AmplitudeC());
    registry.fill(HIST("Events/") + HIST(Stage[stage]) + HIST("hFV0A"), collision.totalFV0AmplitudeA());
    registry.fill(HIST("Events/") + HIST(Stage[stage]) + HIST("hZNEnergy"), znaE, zncE);
    registry.fill(HIST("Events/") + HIST(Stage[stage]) + HIST("hZNTime"), znaT, zncT);
  }

  template <int stage, typename T>
  void fillTrackQA(const T& track)
  {
    const double p = std::sqrt(track.px() * track.px() + track.py() * track.py() + track.pz() * track.pz());
    registry.fill(HIST("Tracks/") + HIST(Stage[stage]) + HIST("hPt"), track.pt());
    registry.fill(HIST("Tracks/") + HIST(Stage[stage]) + HIST("hEta"), eta(track.px(), track.py(), track.pz()));
    registry.fill(HIST("Tracks/") + HIST(Stage[stage]) + HIST("hPhi"), phi(track.px(), track.py()));
    registry.fill(HIST("Tracks/") + HIST(Stage[stage]) + HIST("hDcaXYVsPt"), track.pt(), track.dcaXY());
    registry.fill(HIST("Tracks/") + HIST(Stage[stage]) + HIST("hDcaZ"), track.dcaZ());
    registry.fill(HIST("Tracks/") + HIST(Stage[stage]) + HIST("hItsNCls"), track.itsNCls());
    registry.fill(HIST("Tracks/") + HIST(Stage[stage]) + HIST("hItsChi2NCl"), track.itsChi2NCl());
    registry.fill(HIST("Tracks/") + HIST(Stage[stage]) + HIST("hTpcNClsFound"), track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
    registry.fill(HIST("Tracks/") + HIST(Stage[stage]) + HIST("hTpcCrossedRows"), track.tpcNClsCrossedRows());
    registry.fill(HIST("Tracks/") + HIST(Stage[stage]) + HIST("hTpcChi2NCl"), track.tpcChi2NCl());
    registry.fill(HIST("Tracks/") + HIST(Stage[stage]) + HIST("hTpcSignalVsP"), p, track.tpcSignal());
    registry.fill(HIST("Tracks/") + HIST(Stage[stage]) + HIST("hTpcNSigmaPi"), track.tpcNSigmaPi());
  }

  template <int sel, int chg>
  void fillSystemQA(float m, float pt, float y, float phiSys)
  {
    registry.fill(HIST("System/") + HIST(Selection[sel]) + HIST(Charge[chg]) + HIST("hM"), m);
    registry.fill(HIST("System/") + HIST(Selection[sel]) + HIST(Charge[chg]) + HIST("hPt"), pt);
    registry.fill(HIST("System/") + HIST(Selection[sel]) + HIST(Charge[chg]) + HIST("hPt2"), pt * pt);
    registry.fill(HIST("System/") + HIST(Selection[sel]) + HIST(Charge[chg]) + HIST("hY"), y);
    registry.fill(HIST("System/") + HIST(Selection[sel]) + HIST(Charge[chg]) + HIST("hPhi"), phiSys);
    registry.fill(HIST("System/") + HIST(Selection[sel]) + HIST(Charge[chg]) + HIST("hPtVsM"), m, pt);
  }

  // Contador ACUMULATIVO (igual que el de eventos): cada bin = trazas que pasan hasta ese corte
  template <typename T>
  bool trackPassesCuts(const T& track)
  {
    registry.fill(HIST("Tracks/hSelectionCounter"), 0);
    if (!track.isPVContributor())
      return false;
    registry.fill(HIST("Tracks/hSelectionCounter"), 1);
    if (!track.hasITS())
      return false;
    registry.fill(HIST("Tracks/hSelectionCounter"), 2);
    if (track.itsNCls() < tracksMinItsNClsCut)
      return false;
    registry.fill(HIST("Tracks/hSelectionCounter"), 3);
    if ((track.itsClusterMap() & 0x7) == 0) // al menos un hit en las 3 capas internas
      return false;
    registry.fill(HIST("Tracks/hSelectionCounter"), 4);
    if (track.itsChi2NCl() > itsChi2NClsCut)
      return false;
    registry.fill(HIST("Tracks/hSelectionCounter"), 5);
    if (!track.hasTPC())
      return false;
    registry.fill(HIST("Tracks/hSelectionCounter"), 6);
    if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < tracksMinTpcNClsCut) // clusters ENCONTRADOS
      return false;
    registry.fill(HIST("Tracks/hSelectionCounter"), 7);
    if (track.tpcChi2NCl() < tracksMinTpcChi2NClCut || track.tpcChi2NCl() > tracksMaxTpcChi2NClCut)
      return false;
    registry.fill(HIST("Tracks/hSelectionCounter"), 8);
    if (track.tpcNClsCrossedRows() < tracksMinTpcNClsCrossedRowsCut)
      return false;
    registry.fill(HIST("Tracks/hSelectionCounter"), 9);
    if (track.tpcNClsFindable() <= 0 ||
        static_cast<float>(track.tpcNClsCrossedRows()) / track.tpcNClsFindable() < tracksMinTpcNClsCrossedOverFindableCut)
      return false;
    registry.fill(HIST("Tracks/hSelectionCounter"), 10);
    if (requireTof && !track.hasTOF())
      return false;
    registry.fill(HIST("Tracks/hSelectionCounter"), 11);
    if (track.pt() < tracksMinPtCut)
      return false;
    registry.fill(HIST("Tracks/hSelectionCounter"), 12);
    // DCAxy: parametrizacion estandar de global tracks en O2
    if (std::abs(track.dcaZ()) > dcaZcut || std::abs(track.dcaXY()) > 0.0105 + 0.0350 / std::pow(track.pt(), 1.1))
      return false;
    registry.fill(HIST("Tracks/hSelectionCounter"), 13);
    if (std::abs(eta(track.px(), track.py(), track.pz())) > etaCut)
      return false;
    registry.fill(HIST("Tracks/hSelectionCounter"), 14);
    return true;
  }

  // PID del par: radio en el espacio (nSigma_1, nSigma_2) para cada hipotesis
  template <typename T>
  bool pairPassesPid(const T& selTracks)
  {
    float r2Pi = 0.f, r2El = 0.f, r2Ka = 0.f, r2Pr = 0.f;
    for (const auto& track : selTracks) {
      registry.fill(HIST("Tracks/hTpcNSigmaPiVsEl"), track.tpcNSigmaPi(), track.tpcNSigmaEl());
      r2Pi += track.tpcNSigmaPi() * track.tpcNSigmaPi();
      r2El += track.tpcNSigmaEl() * track.tpcNSigmaEl();
      r2Ka += track.tpcNSigmaKa() * track.tpcNSigmaKa();
      r2Pr += track.tpcNSigmaPr() * track.tpcNSigmaPr();
    }
    registry.fill(HIST("Tracks/hPidRadiusPi"), std::sqrt(r2Pi));
    registry.fill(HIST("Tracks/hPidRadiusEl"), std::sqrt(r2El));

    const float cut = tracksTpcNSigmaPiCut.value;
    if (r2Pi > cut * cut)
      return false;
    if (rejectLowerProbPairs && (r2Pi > r2El || r2Pi > r2Ka || r2Pi > r2Pr))
      return false;
    return true;
  }

  void process(UDCollisions::iterator const& collision, UDtracks const& tracks)
  {
    // ZDC: el SG-producer puede guardar -inf cuando no hay senal
    float znaE = collision.energyCommonZNA(), zncE = collision.energyCommonZNC();
    float znaT = collision.timeZNA(), zncT = collision.timeZNC();
    if (std::isinf(znaE))
      znaE = -999.f;
    if (std::isinf(zncE))
      zncE = -999.f;
    if (std::isinf(znaT))
      znaT = -999.f;
    if (std::isinf(zncT))
      zncT = -999.f;

    fillEventQA<0>(collision, znaE, zncE, znaT, zncT);

    int nPVTracks = 0;
    for (const auto& track : tracks)
      if (track.isPVContributor())
        nPVTracks++;
    registry.fill(HIST("Events/hNumContribVsPVTracks"), nPVTracks, collision.numContrib());

    // ==================== EVENT SELECTION ====================
    registry.fill(HIST("Events/hSelectionCounter"), 0);

    int gap = collision.gapSide();
    if (useTrueGap)
      gap = sgSelector.trueGap(collision, cutTrueGapSideFV0, cutTrueGapSideFT0A, cutTrueGapSideFT0C, cutTrueGapSideZDC);
    if (cutGapSide && gap != gapSide)
      return;
    registry.fill(HIST("Events/hSelectionCounter"), 1);

    if (vtxITSTPCcut && !collision.vtxITSTPC())
      return;
    registry.fill(HIST("Events/hSelectionCounter"), 2);

    if (sbpCut && !collision.sbp())
      return;
    registry.fill(HIST("Events/hSelectionCounter"), 3);

    if (itsROFbCut && !collision.itsROFb())
      return;
    registry.fill(HIST("Events/hSelectionCounter"), 4);

    if (tfbCut && !collision.tfb())
      return;
    registry.fill(HIST("Events/hSelectionCounter"), 5);

    if (useRctFlag && !isGoodRctFlag(collision))
      return;
    registry.fill(HIST("Events/hSelectionCounter"), 6);

    if (std::abs(collision.posZ()) > vZCut)
      return;
    registry.fill(HIST("Events/hSelectionCounter"), 7);

    if (collision.numContrib() != numPVContrib)
      return;
    registry.fill(HIST("Events/hSelectionCounter"), 8);

    if (cutFit && (collision.totalFV0AmplitudeA() > fv0Cut ||
                   collision.totalFT0AmplitudeA() > ft0aCut ||
                   collision.totalFT0AmplitudeC() > ft0cCut))
      return;
    registry.fill(HIST("Events/hSelectionCounter"), 9);

    const int neutronClass = getNeutronClass(znaT, zncT);
    registry.fill(HIST("Events/hNeutronClass"), neutronClass);
    if (neutronClassCut >= 0 && neutronClass != neutronClassCut)
      return;
    registry.fill(HIST("Events/hSelectionCounter"), 10);

    // ==================== TRACK SELECTION ====================
    std::vector<decltype(tracks.begin())> selTracks;
    for (const auto& track : tracks) {
      fillTrackQA<0>(track); // "before" = todas las trazas de eventos que pasaron la seleccion de evento
      if (!trackPassesCuts(track))
        continue;
      selTracks.push_back(track);
    }
    registry.fill(HIST("Tracks/hNSelectedTracks"), selTracks.size());

    // exactamente 2 trazas buenas, ambas contribuyentes del PV y N_contrib == 2
    // -> los dos piones SON el vertice (mismo vertice garantizado)
    if (selTracks.size() != 2)
      return;
    registry.fill(HIST("Events/hSelectionCounter"), 11);

    if (!pairPassesPid(selTracks))
      return;
    registry.fill(HIST("Events/hSelectionCounter"), 12);

    for (const auto& track : selTracks)
      fillTrackQA<1>(track);
    fillEventQA<1>(collision, znaE, zncE, znaT, zncT);

    // ==================== SYSTEM ====================
    const auto& t0 = selTracks[0];
    const auto& t1 = selTracks[1];
    ROOT::Math::PxPyPzMVector p0(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassPionCharged);
    ROOT::Math::PxPyPzMVector p1(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
    ROOT::Math::PxPyPzMVector system = p0 + p1;

    const int totalCharge = t0.sign() + t1.sign();
    const float m = system.M();
    const float pt = system.Pt();
    const float y = system.Rapidity();
    const float phiSys = system.Phi() + o2::constants::math::PI; // [0, 2pi)

    const bool passM = m >= systemMassMinCut && m <= systemMassMaxCut;
    const bool passPt = pt <= systemPtCut;
    const bool passY = std::abs(y) <= systemYCut;
    const bool selected = passM && passPt && passY;

    if (totalCharge == 0) {
      registry.fill(HIST("Events/hSelectionCounter"), 13);
      fillSystemQA<0, 0>(m, pt, y, phiSys);
      if (passPt && passY)
        registry.fill(HIST("System/nMinus1/hM"), m);
      if (passM && passY) {
        registry.fill(HIST("System/nMinus1/hPt"), pt);
        registry.fill(HIST("System/nMinus1/hPt2"), pt * pt);
      }
      if (selected) {
        registry.fill(HIST("Events/hSelectionCounter"), 14);
        fillSystemQA<1, 0>(m, pt, y, phiSys);
      }
    } else {
      // like-sign (++ o --): control de fondo combinatorio, NO va al arbol
      fillSystemQA<0, 1>(m, pt, y, phiSys);
      if (selected)
        fillSystemQA<1, 1>(m, pt, y, phiSys);
      return;
    }

    // ==================== FILL TREE (solo pi+ pi- del mismo vertice) ====================
    if (saveOnlySelected) {
      if (!selected)
        return;
    } else {
      if (!passY || pt > treeMaxPtCut)
        return;
    }

    auto piPlus = (t0.sign() > 0) ? t0 : t1;
    auto piMinus = (t0.sign() > 0) ? t1 : t0;

    systemTree(collision.runNumber(), collision.flags(),
               collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(),
               collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), collision.totalFV0AmplitudeA(),
               znaE, zncE, znaT, zncT, neutronClass,
               m, pt, y, phiSys, selected,
               piPlus.pt(), eta(piPlus.px(), piPlus.py(), piPlus.pz()), phi(piPlus.px(), piPlus.py()),
               piPlus.tpcNSigmaPi(), piPlus.tpcNSigmaEl(),
               piMinus.pt(), eta(piMinus.px(), piMinus.py(), piMinus.pz()), phi(piMinus.px(), piMinus.py()),
               piMinus.tpcNSigmaPi(), piMinus.tpcNSigmaEl());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcTwoPionAnalysis>(cfgc, TaskName{"upc-two-pion-analysis"})};
}
