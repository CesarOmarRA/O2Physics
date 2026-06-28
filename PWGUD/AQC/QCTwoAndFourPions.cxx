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
/// \brief  Quality Control task for 2 and 4 pions using numContrib/2 on AO2D standard tables

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>

#include <cmath>
#include <vector>

using namespace o2;
using namespace o2::framework;

// ============ DEFINICIÓN DE TABLAS ESTÁNDAR AO2D ============
using Collisions = soa::Join<aod::Collisions, aod::EvSels>;
using Tracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, 
                         aod::TracksDCA, aod::TrackSelectionExtension, 
                         aod::pidTPCFullPi>;

struct QCTwoAndFourPions {
  
  // ============ CONFIGURABLES (CORTES DE TU CÓDIGO ORIGINAL) ============
  
  // Cortes del sistema
  Configurable<double> systemYCut{"systemYCut", 0.5, "Max Rapidity of system"};
  Configurable<double> systemPtCut{"systemPtCut", 0.1, "Min Pt of system"};
  
  // Cortes de masa para 2 piones
  Configurable<double> system2PiMassMinCut{"system2PiMassMinCut", 0.5, "Min Mass for 2#pi system"};
  Configurable<double> system2PiMassMaxCut{"system2PiMassMaxCut", 1.0, "Max Mass for 2#pi system"};
  
  // Cortes de masa para 4 piones
  Configurable<double> system4PiMassMinCut{"system4PiMassMinCut", 0.8, "Min Mass for 4#pi system"};
  Configurable<double> system4PiMassMaxCut{"system4PiMassMaxCut", 2.2, "Max Mass for 4#pi system"};
  
  // Cortes de tracks
  Configurable<double> etaCut{"etaCut", 0.9, "Track Pseudorapidity"};
  Configurable<float> minPtTrack{"minPtTrack", 0.1, "Minimum track pT"};
  
  // Cortes de evento
  Configurable<float> vZCut{"vZCut", 10.0, "Cut on vertex Z position"};
  
  // Cortes de calidad de tracks
  Configurable<bool> useOnlyPVtracks{"useOnlyPVtracks", true, "Use only PV tracks"};
  Configurable<float> tpcChi2NClsCut{"tpcChi2NClsCut", 5.0, "TPC chi2/N clusters cut"};
  Configurable<float> nSigmaTPCcut{"nSigmaTPCcut", 5.0, "TPC nSigma cut for pions"};
  Configurable<float> dcaZcut{"dcaZcut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXYcut{"dcaXYcut", 0.024, "dcaXY cut"};
  Configurable<int> minTPCFindableClusters{"minTPCFindableClusters", 70, "Minimum TPC findable clusters"};

  // ============ HISTOGRAMAS ============
  HistogramRegistry registry{
    "registry",
    {
        // Event flow
        {"Events/Flow", "Event flow;Cut;Counts", {HistType::kTH1F, {{10, 0, 10}}}},
        {"Events/VertexZ", "Vertex Z;z (cm);Counts", {HistType::kTH1F, {{200, -20, 20}}}},
        {"Events/NumContrib", "Number of contributors;N_{contrib};Counts", {HistType::kTH1F, {{50, 0, 50}}}},
        
        // Track quality
        {"Tracks/Pt", "Track p_{T};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{200, 0, 2}}}},
        {"Tracks/Eta", "Track #eta;#eta;Counts", {HistType::kTH1F, {{200, -2, 2}}}},
        {"Tracks/TPCNSigmaPi", "TPC n#sigma for #pi;n#sigma;Counts", {HistType::kTH1F, {{200, -10, 10}}}},
        {"Tracks/Charge", "Track charge;Charge;Counts", {HistType::kTH1F, {{3, -1.5, 1.5}}}},
        
        // Sistema de 2 piones
        {"System2Pi/hM", "2#pi Invariant Mass;m (GeV/#it{c}^{2});Counts", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
        {"System2Pi/hPt", "2#pi p_{T};p_{T} (GeV/#it{c});Counts", {HistType::kTH1F, {{1000, 0.0, 2.0}}}},
        {"System2Pi/hEta", "2#pi #eta;#eta;Counts", {HistType::kTH1F, {{180, -0.9, 0.9}}}},
        {"System2Pi/hPhi", "2#pi #phi;#phi;Counts", {HistType::kTH1F, {{180, 0.0, 6.28}}}},
        {"System2Pi/hY", "2#pi rapidity;y;Counts", {HistType::kTH1F, {{180, -0.9, 0.9}}}},
        {"System2Pi/hMassVsPt", "2#pi Mass vs p_{T};m (GeV/c^{2});p_{T} (GeV/c)", {HistType::kTH2F, {{200, 0, 5}, {200, 0, 2}}}},
        
        // Sistema de 4 piones
        {"System4Pi/hM", "4#pi Invariant Mass;m (GeV/#it{c}^{2});Counts", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
        {"System4Pi/hPt", "4#pi p_{T};p_{T} (GeV/#it{c});Counts", {HistType::kTH1F, {{1000, 0.0, 2.0}}}},
        {"System4Pi/hEta", "4#pi #eta;#eta;Counts", {HistType::kTH1F, {{180, -0.9, 0.9}}}},
        {"System4Pi/hPhi", "4#pi #phi;#phi;Counts", {HistType::kTH1F, {{180, 0.0, 6.28}}}},
        {"System4Pi/hY", "4#pi rapidity;y;Counts", {HistType::kTH1F, {{180, -0.9, 0.9}}}},
        {"System4Pi/hMassVsPt", "4#pi Mass vs p_{T};m (GeV/c^{2});p_{T} (GeV/c)", {HistType::kTH2F, {{200, 0, 5}, {200, 0, 2}}}},
        
        // Comparación
        {"Comparison/MassBoth", "Mass comparison;m (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0, 5}}}},
        {"Comparison/PtBoth", "p_{T} comparison;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0, 2}}}},
        
        // QC Monitoring
        {"QC/NumContribDistribution", "numContrib distribution;numContrib;Counts", {HistType::kTH1F, {{50, 0, 50}}}},
        {"QC/PionsFound", "Pions found per event;N_{#pi};Counts", {HistType::kTH1F, {{10, 0, 10}}}},
    }};

  void init(InitContext&)
  {
    auto hFlow = registry.get<TH1>(HIST("Events/Flow"));
    hFlow->GetXaxis()->SetBinLabel(1, "All events");
    hFlow->GetXaxis()->SetBinLabel(2, "Vertex Z cut");
    hFlow->GetXaxis()->SetBinLabel(3, "Has tracks");
    hFlow->GetXaxis()->SetBinLabel(4, "Found pions");
    hFlow->GetXaxis()->SetBinLabel(5, "Mass cut passed");
  }

  void process(Collisions::iterator const& collision, Tracks const& tracks)
  {
    // ============ FLOW ============
    registry.fill(HIST("Events/Flow"), 0);
    
    // ============ SELECCIÓN DE EVENTOS ============
    // Vertex Z cut
    if (std::abs(collision.posZ()) > vZCut) return;
    registry.fill(HIST("Events/Flow"), 1);
    
    // Que haya tracks
    if (tracks.size() == 0) return;
    registry.fill(HIST("Events/Flow"), 2);
    
    // Llenar histogramas de evento
    registry.fill(HIST("Events/VertexZ"), collision.posZ());
    registry.fill(HIST("Events/NumContrib"), collision.numContrib());
    registry.fill(HIST("QC/NumContribDistribution"), collision.numContrib());
    
    // ============ DETERMINAR NÚMERO DE PIONES USANDO numContrib/2 ============
    int numContrib = collision.numContrib();
    int expectedPairs = numContrib / 2;
    int expectedPions = expectedPairs * 2;
    
    // Solo procesar si es 2 o 4 piones
    if (expectedPions != 2 && expectedPions != 4) return;
    
    // ============ SELECCIÓN DE TRACKS ============
    std::vector<decltype(tracks.begin())> posPions;
    std::vector<decltype(tracks.begin())> negPions;
    posPions.reserve(expectedPairs);
    negPions.reserve(expectedPairs);
    
    for (const auto& track : tracks) {
      // PV contributor
      if (useOnlyPVtracks && !track.isPVContributor()) continue;
      
      // Calidad básica
      if (!track.hasITS() || !track.hasTPC()) continue;
      
      // Cortes cinemáticos
      if (track.pt() < minPtTrack) continue;
      if (std::abs(track.eta()) > etaCut) continue;
      
      // Calidad TPC
      if (track.tpcChi2NCl() > tpcChi2NClsCut) continue;
      if (track.tpcNClsFindable() < minTPCFindableClusters) continue;
      
      // PID: nSigma para piones
      if (std::abs(track.tpcNSigmaPi()) > nSigmaTPCcut) continue;
      
      // DCA cuts
      if (std::abs(track.dcaZ()) > dcaZcut) continue;
      
      float maxDCAxy = 0.0105 + 0.035 / std::pow(track.pt(), 1.1);
      if (dcaXYcut == 0 && std::abs(track.dcaXY()) > maxDCAxy) continue;
      else if (dcaXYcut > 0 && std::abs(track.dcaXY()) > dcaXYcut) continue;
      
      // Llenar histogramas de tracks
      registry.fill(HIST("Tracks/Pt"), track.pt());
      registry.fill(HIST("Tracks/Eta"), track.eta());
      registry.fill(HIST("Tracks/TPCNSigmaPi"), track.tpcNSigmaPi());
      registry.fill(HIST("Tracks/Charge"), track.sign());
      
      // Seleccionar piones según carga
      if (track.sign() > 0 && posPions.size() < expectedPairs) {
        posPions.push_back(track);
      } else if (track.sign() < 0 && negPions.size() < expectedPairs) {
        negPions.push_back(track);
      }
      
      if (posPions.size() == expectedPairs && negPions.size() == expectedPairs) break;
    }
    
    // Verificar que encontramos los piones esperados
    if (posPions.size() != expectedPairs || negPions.size() != expectedPairs) return;
    
    registry.fill(HIST("QC/PionsFound"), expectedPions);
    registry.fill(HIST("Events/Flow"), 3);
    
    // ============ RECONSTRUCCIÓN DEL SISTEMA ============
    std::vector<decltype(tracks.begin())> selectedTracks;
    selectedTracks.insert(selectedTracks.end(), posPions.begin(), posPions.end());
    selectedTracks.insert(selectedTracks.end(), negPions.begin(), negPions.end());
    
    TLorentzVector system;
    for (const auto& track : selectedTracks) {
      TLorentzVector pionVec;
      pionVec.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), 
                           o2::constants::physics::MassPionCharged);
      system += pionVec;
    }
    
    // ============ APLICAR CORTES DE MASA SEGÚN NÚMERO DE PIONES ============
    bool passMassCut = false;
    
    if (expectedPions == 2) {
      if (system.M() >= system2PiMassMinCut && system.M() <= system2PiMassMaxCut) {
        passMassCut = true;
      }
    } else if (expectedPions == 4) {
      if (system.M() >= system4PiMassMinCut && system.M() <= system4PiMassMaxCut) {
        passMassCut = true;
      }
    }
    
    if (!passMassCut) return;
    registry.fill(HIST("Events/Flow"), 4);
    
    // Cortes comunes de Pt y Rapidez
    if (system.Pt() < systemPtCut) return;
    if (std::abs(system.Rapidity()) > systemYCut) return;
    
    // ============ LLENAR HISTOGRAMAS ============
    if (expectedPions == 2) {
      registry.fill(HIST("System2Pi/hM"), system.M());
      registry.fill(HIST("System2Pi/hPt"), system.Pt());
      registry.fill(HIST("System2Pi/hEta"), system.Eta());
      registry.fill(HIST("System2Pi/hPhi"), system.Phi());
      registry.fill(HIST("System2Pi/hY"), system.Rapidity());
      registry.fill(HIST("System2Pi/hMassVsPt"), system.M(), system.Pt());
      registry.fill(HIST("Comparison/MassBoth"), system.M());
      registry.fill(HIST("Comparison/PtBoth"), system.Pt());
    } 
    else if (expectedPions == 4) {
      registry.fill(HIST("System4Pi/hM"), system.M());
      registry.fill(HIST("System4Pi/hPt"), system.Pt());
      registry.fill(HIST("System4Pi/hEta"), system.Eta());
      registry.fill(HIST("System4Pi/hPhi"), system.Phi());
      registry.fill(HIST("System4Pi/hY"), system.Rapidity());
      registry.fill(HIST("System4Pi/hMassVsPt"), system.M(), system.Pt());
      registry.fill(HIST("Comparison/MassBoth"), system.M());
      registry.fill(HIST("Comparison/PtBoth"), system.Pt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<QCTwoAndFourPions>(cfgc)};
}