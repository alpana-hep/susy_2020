//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 17 04:59:05 2024 by ROOT version 6.18/04
// from TTree PreSelection/PreSelection
// found on file: root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2018A-UL2018-v3/SingleMuon/2585_RA2AnalysisTree.root
//////////////////////////////////////////////////////////

#ifndef allin1_h
#define allin1_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "Math/GenVector/PtEtaPhiE4D.h"
#include "vector"
#include "vector"
#include "vector"

class allin1 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxElectrons = 3;
   static constexpr Int_t kMaxHLTMuonObjects = 2;
   static constexpr Int_t kMaxJets = 36;
   static constexpr Int_t kMaxJetsAK15 = 9;
   static constexpr Int_t kMaxJetsAK15_subjets = 10;
   static constexpr Int_t kMaxJetsAK8 = 21;
   static constexpr Int_t kMaxJetsAK8_subjets = 32;
   static constexpr Int_t kMaxJetsConstituents = 386;
   static constexpr Int_t kMaxMuons = 4;
   static constexpr Int_t kMaxPhotons = 2;
   static constexpr Int_t kMaxTAPElectronTracks = 4;
   static constexpr Int_t kMaxTAPMuonTracks = 4;
   static constexpr Int_t kMaxTAPPionTracks = 27;

   // Declaration of leaf types
   UInt_t          RunNum;
   UInt_t          LumiBlockNum;
   ULong64_t       EvtNum;
   Int_t           BadChargedCandidateFilter;
   Bool_t          BadPFMuonDzFilter;
   Int_t           BadPFMuonFilter;
   Int_t           BTagsDeepCSV;
   Float_t         CaloMET;
   Float_t         CaloMETPhi;
   Int_t           CSCTightHaloFilter;
   Float_t         DeltaPhi1;
   Float_t         DeltaPhi1_AK8;
   Float_t         DeltaPhi2;
   Float_t         DeltaPhi2_AK8;
   Float_t         DeltaPhi3;
   Float_t         DeltaPhi4;
   Float_t         DeltaPhiMin_AK8;
   Int_t           ecalBadCalibFilter;
   Int_t           EcalDeadCellBoundaryEnergyFilter;
   Int_t           EcalDeadCellTriggerPrimitiveFilter;
   Int_t           eeBadScFilter;
   Int_t           Electrons_;
   Float_t         Electrons_fCoordinates_fPt[kMaxElectrons];   //[Electrons_]
   Float_t         Electrons_fCoordinates_fEta[kMaxElectrons];   //[Electrons_]
   Float_t         Electrons_fCoordinates_fPhi[kMaxElectrons];   //[Electrons_]
   Float_t         Electrons_fCoordinates_fE[kMaxElectrons];   //[Electrons_]
   vector<int>     *Electrons_charge;
   vector<float>   *Electrons_iso;
   vector<bool>    *Electrons_mediumID;
   vector<float>   *Electrons_MTW;
   vector<bool>    *Electrons_passIso;
   vector<bool>    *Electrons_tightID;
   Float_t         fixedGridRhoFastjetAll;
   Float_t         GenMT2_AK8;
   Int_t           globalSuperTightHalo2016Filter;
   Int_t           globalTightHalo2016Filter;
   Bool_t          hasGenPromptPhoton;
   Int_t           HBHEIsoNoiseFilter;
   Int_t           HBHENoiseFilter;
   Int_t           hfNoisyHitsFilter;
   Int_t           HLTMuonObjects_;
   Float_t         HLTMuonObjects_fCoordinates_fPt[kMaxHLTMuonObjects];   //[HLTMuonObjects_]
   Float_t         HLTMuonObjects_fCoordinates_fEta[kMaxHLTMuonObjects];   //[HLTMuonObjects_]
   Float_t         HLTMuonObjects_fCoordinates_fPhi[kMaxHLTMuonObjects];   //[HLTMuonObjects_]
   Float_t         HLTMuonObjects_fCoordinates_fE[kMaxHLTMuonObjects];   //[HLTMuonObjects_]
   Float_t         HT;
   Float_t         HT5;
   Float_t         HTOnline;
   Int_t           isoElectronTracks;
   Int_t           isoMuonTracks;
   Int_t           isoPionTracks;
   Bool_t          JetID;
   Bool_t          JetIDAK15;
   Bool_t          JetIDAK8;
   Int_t           Jets_;
   Float_t         Jets_fCoordinates_fPt[kMaxJets];   //[Jets_]
   Float_t         Jets_fCoordinates_fEta[kMaxJets];   //[Jets_]
   Float_t         Jets_fCoordinates_fPhi[kMaxJets];   //[Jets_]
   Float_t         Jets_fCoordinates_fE[kMaxJets];   //[Jets_]
   vector<float>   *Jets_axismajor;
   vector<float>   *Jets_axisminor;
   vector<float>   *Jets_bDiscriminatorCSV;
   vector<float>   *Jets_bJetTagDeepCSVBvsAll;
   vector<float>   *Jets_bJetTagDeepCSVprobb;
   vector<float>   *Jets_bJetTagDeepCSVprobbb;
   vector<float>   *Jets_bJetTagDeepCSVprobc;
   vector<float>   *Jets_bJetTagDeepCSVprobudsg;
   vector<float>   *Jets_bJetTagDeepFlavourprobb;
   vector<float>   *Jets_bJetTagDeepFlavourprobbb;
   vector<float>   *Jets_bJetTagDeepFlavourprobc;
   vector<float>   *Jets_bJetTagDeepFlavourprobg;
   vector<float>   *Jets_bJetTagDeepFlavourproblepb;
   vector<float>   *Jets_bJetTagDeepFlavourprobuds;
   vector<float>   *Jets_chargedEmEnergyFraction;
   vector<float>   *Jets_chargedHadronEnergyFraction;
   vector<int>     *Jets_chargedHadronMultiplicity;
   vector<int>     *Jets_chargedMultiplicity;
   vector<float>   *Jets_electronEnergyFraction;
   vector<int>     *Jets_electronMultiplicity;
   vector<int>     *Jets_hadronFlavor;
   vector<float>   *Jets_hfEMEnergyFraction;
   vector<float>   *Jets_hfHadronEnergyFraction;
   vector<bool>    *Jets_HTMask;
   vector<bool>    *Jets_ID;
   vector<float>   *Jets_jecFactor;
   vector<float>   *Jets_jecUnc;
   vector<bool>    *Jets_LeptonMask;
   vector<bool>    *Jets_MHTMask;
   vector<int>     *Jets_multiplicity;
   vector<float>   *Jets_muonEnergyFraction;
   vector<int>     *Jets_muonMultiplicity;
   vector<float>   *Jets_neutralEmEnergyFraction;
   vector<float>   *Jets_neutralHadronEnergyFraction;
   vector<int>     *Jets_neutralHadronMultiplicity;
   vector<int>     *Jets_neutralMultiplicity;
   vector<int>     *Jets_partonFlavor;
   vector<float>   *Jets_photonEnergyFraction;
   vector<int>     *Jets_photonMultiplicity;
   vector<float>   *Jets_pileupId;
   vector<float>   *Jets_ptD;
   vector<float>   *Jets_qgLikelihood;
   Int_t           JetsAK15_;
   Float_t         JetsAK15_fCoordinates_fPt[kMaxJetsAK15];   //[JetsAK15_]
   Float_t         JetsAK15_fCoordinates_fEta[kMaxJetsAK15];   //[JetsAK15_]
   Float_t         JetsAK15_fCoordinates_fPhi[kMaxJetsAK15];   //[JetsAK15_]
   Float_t         JetsAK15_fCoordinates_fE[kMaxJetsAK15];   //[JetsAK15_]
   vector<float>   *JetsAK15_axismajor;
   vector<float>   *JetsAK15_axisminor;
   vector<float>   *JetsAK15_chargedEmEnergyFraction;
   vector<float>   *JetsAK15_chargedHadronEnergyFraction;
   vector<int>     *JetsAK15_chargedHadronMultiplicity;
   vector<int>     *JetsAK15_chargedMultiplicity;
   vector<int>     *JetsAK15_constituentsIndex;
   vector<int>     *JetsAK15_constituentsIndexCounts;
   vector<float>   *JetsAK15_DeepMassDecorrelTagbbvsLight;
   vector<float>   *JetsAK15_DeepMassDecorrelTagHbbvsQCD;
   vector<float>   *JetsAK15_DeepMassDecorrelTagTvsQCD;
   vector<float>   *JetsAK15_DeepMassDecorrelTagWvsQCD;
   vector<float>   *JetsAK15_DeepMassDecorrelTagZbbvsQCD;
   vector<float>   *JetsAK15_DeepMassDecorrelTagZHbbvsQCD;
   vector<float>   *JetsAK15_DeepMassDecorrelTagZvsQCD;
   vector<float>   *JetsAK15_DeepTagHbbvsQCD;
   vector<float>   *JetsAK15_DeepTagTvsQCD;
   vector<float>   *JetsAK15_DeepTagWvsQCD;
   vector<float>   *JetsAK15_DeepTagZbbvsQCD;
   vector<float>   *JetsAK15_DeepTagZvsQCD;
   vector<float>   *JetsAK15_doubleBDiscriminator;
   vector<float>   *JetsAK15_ecfC2b1;
   vector<float>   *JetsAK15_ecfC2b2;
   vector<float>   *JetsAK15_ecfD2b1;
   vector<float>   *JetsAK15_ecfD2b2;
   vector<float>   *JetsAK15_ecfM2b1;
   vector<float>   *JetsAK15_ecfM2b2;
   vector<float>   *JetsAK15_ecfN2b1;
   vector<float>   *JetsAK15_ecfN2b2;
   vector<float>   *JetsAK15_electronEnergyFraction;
   vector<int>     *JetsAK15_electronMultiplicity;
   vector<float>   *JetsAK15_girth;
   vector<int>     *JetsAK15_hadronFlavor;
   vector<float>   *JetsAK15_hfEMEnergyFraction;
   vector<float>   *JetsAK15_hfHadronEnergyFraction;
   vector<bool>    *JetsAK15_ID;
   vector<float>   *JetsAK15_jecFactor;
   vector<int>     *JetsAK15_multiplicity;
   vector<float>   *JetsAK15_muonEnergyFraction;
   vector<int>     *JetsAK15_muonMultiplicity;
   vector<float>   *JetsAK15_neutralEmEnergyFraction;
   vector<float>   *JetsAK15_neutralHadronEnergyFraction;
   vector<float>   *JetsAK15_neutralHadronMultiplicity;
   vector<float>   *JetsAK15_neutralMultiplicity;
   vector<float>   *JetsAK15_NsubjettinessTau1;
   vector<float>   *JetsAK15_NsubjettinessTau2;
   vector<float>   *JetsAK15_NsubjettinessTau3;
   vector<float>   *JetsAK15_NsubjettinessTau4;
   vector<int>     *JetsAK15_NumBhadrons;
   vector<int>     *JetsAK15_NumChadrons;
   vector<int>     *JetsAK15_partonFlavor;
   vector<float>   *JetsAK15_pfMassIndependentDeepDoubleBvLJetTagsProbHbb;
   vector<float>   *JetsAK15_photonEnergyFraction;
   vector<float>   *JetsAK15_photonMultiplicity;
   vector<float>   *JetsAK15_ptD;
   vector<float>   *JetsAK15_softDropMass;
   vector<float>   *JetsAK15_softDropMassBeta1;
   Int_t           JetsAK15_subjets_;
   Float_t         JetsAK15_subjets_fCoordinates_fPt[kMaxJetsAK15_subjets];   //[JetsAK15_subjets_]
   Float_t         JetsAK15_subjets_fCoordinates_fEta[kMaxJetsAK15_subjets];   //[JetsAK15_subjets_]
   Float_t         JetsAK15_subjets_fCoordinates_fPhi[kMaxJetsAK15_subjets];   //[JetsAK15_subjets_]
   Float_t         JetsAK15_subjets_fCoordinates_fE[kMaxJetsAK15_subjets];   //[JetsAK15_subjets_]
   vector<int>     *JetsAK15_subjetsCounts;
   Int_t           JetsAK8_;
   Float_t         JetsAK8_fCoordinates_fPt[kMaxJetsAK8];   //[JetsAK8_]
   Float_t         JetsAK8_fCoordinates_fEta[kMaxJetsAK8];   //[JetsAK8_]
   Float_t         JetsAK8_fCoordinates_fPhi[kMaxJetsAK8];   //[JetsAK8_]
   Float_t         JetsAK8_fCoordinates_fE[kMaxJetsAK8];   //[JetsAK8_]
   vector<float>   *JetsAK8_axismajor;
   vector<float>   *JetsAK8_axisminor;
   vector<float>   *JetsAK8_chargedEmEnergyFraction;
   vector<float>   *JetsAK8_chargedHadronEnergyFraction;
   vector<int>     *JetsAK8_chargedHadronMultiplicity;
   vector<int>     *JetsAK8_chargedMultiplicity;
   vector<int>     *JetsAK8_constituentsIndex;
   vector<int>     *JetsAK8_constituentsIndexCounts;
   vector<float>   *JetsAK8_DeepMassDecorrelTagbbvsLight;
   vector<float>   *JetsAK8_DeepMassDecorrelTagHbbvsQCD;
   vector<float>   *JetsAK8_DeepMassDecorrelTagTvsQCD;
   vector<float>   *JetsAK8_DeepMassDecorrelTagWvsQCD;
   vector<float>   *JetsAK8_DeepMassDecorrelTagZbbvsQCD;
   vector<float>   *JetsAK8_DeepMassDecorrelTagZHbbvsQCD;
   vector<float>   *JetsAK8_DeepMassDecorrelTagZvsQCD;
   vector<float>   *JetsAK8_DeepTagHbbvsQCD;
   vector<float>   *JetsAK8_DeepTagTvsQCD;
   vector<float>   *JetsAK8_DeepTagWvsQCD;
   vector<float>   *JetsAK8_DeepTagZbbvsQCD;
   vector<float>   *JetsAK8_DeepTagZvsQCD;
   vector<float>   *JetsAK8_doubleBDiscriminator;
   vector<float>   *JetsAK8_ecfN2b1;
   vector<float>   *JetsAK8_ecfN2b2;
   vector<float>   *JetsAK8_ecfN3b1;
   vector<float>   *JetsAK8_ecfN3b2;
   vector<float>   *JetsAK8_electronEnergyFraction;
   vector<int>     *JetsAK8_electronMultiplicity;
   vector<float>   *JetsAK8_girth;
   vector<int>     *JetsAK8_hadronFlavor;
   vector<float>   *JetsAK8_hfEMEnergyFraction;
   vector<float>   *JetsAK8_hfHadronEnergyFraction;
   vector<bool>    *JetsAK8_ID;
   vector<float>   *JetsAK8_jecFactor;
   vector<float>   *JetsAK8_jecUnc;
   vector<int>     *JetsAK8_multiplicity;
   vector<float>   *JetsAK8_muonEnergyFraction;
   vector<int>     *JetsAK8_muonMultiplicity;
   vector<float>   *JetsAK8_neutralEmEnergyFraction;
   vector<float>   *JetsAK8_neutralHadronEnergyFraction;
   vector<float>   *JetsAK8_neutralHadronMultiplicity;
   vector<float>   *JetsAK8_neutralMultiplicity;
   vector<float>   *JetsAK8_NsubjettinessTau1;
   vector<float>   *JetsAK8_NsubjettinessTau2;
   vector<float>   *JetsAK8_NsubjettinessTau3;
   vector<int>     *JetsAK8_NumBhadrons;
   vector<int>     *JetsAK8_NumChadrons;
   vector<int>     *JetsAK8_partonFlavor;
   vector<float>   *JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb;
   vector<float>   *JetsAK8_photonEnergyFraction;
   vector<float>   *JetsAK8_photonMultiplicity;
   vector<float>   *JetsAK8_ptD;
   vector<float>   *JetsAK8_softDropMass;
   Int_t           JetsAK8_subjets_;
   Float_t         JetsAK8_subjets_fCoordinates_fPt[kMaxJetsAK8_subjets];   //[JetsAK8_subjets_]
   Float_t         JetsAK8_subjets_fCoordinates_fEta[kMaxJetsAK8_subjets];   //[JetsAK8_subjets_]
   Float_t         JetsAK8_subjets_fCoordinates_fPhi[kMaxJetsAK8_subjets];   //[JetsAK8_subjets_]
   Float_t         JetsAK8_subjets_fCoordinates_fE[kMaxJetsAK8_subjets];   //[JetsAK8_subjets_]
   vector<int>     *JetsAK8_subjetsCounts;
   vector<float>   *JetsAK8_subjets_axismajor;
   vector<float>   *JetsAK8_subjets_axisminor;
   vector<float>   *JetsAK8_subjets_jecFactor;
   vector<int>     *JetsAK8_subjets_multiplicity;
   vector<float>   *JetsAK8_subjets_ptD;
   Int_t           JetsConstituents_;
   Float_t         JetsConstituents_fCoordinates_fPt[kMaxJetsConstituents];   //[JetsConstituents_]
   Float_t         JetsConstituents_fCoordinates_fEta[kMaxJetsConstituents];   //[JetsConstituents_]
   Float_t         JetsConstituents_fCoordinates_fPhi[kMaxJetsConstituents];   //[JetsConstituents_]
   Float_t         JetsConstituents_fCoordinates_fE[kMaxJetsConstituents];   //[JetsConstituents_]
   vector<float>   *JetsConstituents_dxy;
   vector<float>   *JetsConstituents_dxysig;
   vector<float>   *JetsConstituents_dz;
   vector<float>   *JetsConstituents_dzsig;
   vector<int>     *JetsConstituents_PdgId;
   vector<float>   *JetsConstituents_PuppiWeight;
   Float_t         MET;
   Float_t         METPhi;
   Float_t         METSignificance;
   Float_t         MHT;
   Float_t         MHTOnline;
   Float_t         MHTPhi;
   Float_t         MJJ_AK8;
   Float_t         Mmc_AK8;
   Float_t         MT_AK8;
   Int_t           Muons_;
   Float_t         Muons_fCoordinates_fPt[kMaxMuons];   //[Muons_]
   Float_t         Muons_fCoordinates_fEta[kMaxMuons];   //[Muons_]
   Float_t         Muons_fCoordinates_fPhi[kMaxMuons];   //[Muons_]
   Float_t         Muons_fCoordinates_fE[kMaxMuons];   //[Muons_]
   vector<int>     *Muons_charge;
   vector<float>   *Muons_iso;
   vector<bool>    *Muons_mediumID;
   vector<float>   *Muons_MTW;
   vector<bool>    *Muons_passIso;
   vector<bool>    *Muons_tightID;
   Int_t           nAllVertices;
   Int_t           NElectrons;
   Int_t           NJets;
   Int_t           NMuons;
   Float_t         NonPrefiringProb;
   Float_t         NonPrefiringProbDown;
   Float_t         NonPrefiringProbECAL;
   Float_t         NonPrefiringProbECALDown;
   Float_t         NonPrefiringProbECALUp;
   Float_t         NonPrefiringProbMuon;
   Float_t         NonPrefiringProbMuonDown;
   Float_t         NonPrefiringProbMuonUp;
   Float_t         NonPrefiringProbUp;
   Int_t           NVtx;
   Float_t         PFCaloMETRatio;
   Int_t           Photons_;
   Float_t         Photons_fCoordinates_fPt[kMaxPhotons];   //[Photons_]
   Float_t         Photons_fCoordinates_fEta[kMaxPhotons];   //[Photons_]
   Float_t         Photons_fCoordinates_fPhi[kMaxPhotons];   //[Photons_]
   Float_t         Photons_fCoordinates_fE[kMaxPhotons];   //[Photons_]
   vector<bool>    *Photons_electronFakes;
   vector<bool>    *Photons_fullID;
   vector<float>   *Photons_genMatched;
   vector<float>   *Photons_hadTowOverEM;
   vector<bool>    *Photons_hasPixelSeed;
   vector<float>   *Photons_isEB;
   vector<bool>    *Photons_nonPrompt;
   vector<float>   *Photons_passElectronVeto;
   vector<float>   *Photons_pfChargedIso;
   vector<float>   *Photons_pfChargedIsoRhoCorr;
   vector<float>   *Photons_pfGammaIso;
   vector<float>   *Photons_pfGammaIsoRhoCorr;
   vector<float>   *Photons_pfNeutralIso;
   vector<float>   *Photons_pfNeutralIsoRhoCorr;
   vector<float>   *Photons_sigmaIetaIeta;
   Float_t         PrescaleWeightHT;
   Int_t           PrimaryVertexFilter;
   Int_t           TAPElectronTracks_;
   Float_t         TAPElectronTracks_fCoordinates_fPt[kMaxTAPElectronTracks];   //[TAPElectronTracks_]
   Float_t         TAPElectronTracks_fCoordinates_fEta[kMaxTAPElectronTracks];   //[TAPElectronTracks_]
   Float_t         TAPElectronTracks_fCoordinates_fPhi[kMaxTAPElectronTracks];   //[TAPElectronTracks_]
   Float_t         TAPElectronTracks_fCoordinates_fE[kMaxTAPElectronTracks];   //[TAPElectronTracks_]
   vector<float>   *TAPElectronTracks_dxypv;
   vector<bool>    *TAPElectronTracks_leptonMatch;
   vector<float>   *TAPElectronTracks_mT;
   vector<float>   *TAPElectronTracks_pfRelIso03chg;
   vector<float>   *TAPElectronTracks_trkiso;
   Int_t           TAPMuonTracks_;
   Float_t         TAPMuonTracks_fCoordinates_fPt[kMaxTAPMuonTracks];   //[TAPMuonTracks_]
   Float_t         TAPMuonTracks_fCoordinates_fEta[kMaxTAPMuonTracks];   //[TAPMuonTracks_]
   Float_t         TAPMuonTracks_fCoordinates_fPhi[kMaxTAPMuonTracks];   //[TAPMuonTracks_]
   Float_t         TAPMuonTracks_fCoordinates_fE[kMaxTAPMuonTracks];   //[TAPMuonTracks_]
   vector<float>   *TAPMuonTracks_dxypv;
   vector<bool>    *TAPMuonTracks_leptonMatch;
   vector<float>   *TAPMuonTracks_mT;
   vector<float>   *TAPMuonTracks_pfRelIso03chg;
   vector<float>   *TAPMuonTracks_trkiso;
   Int_t           TAPPionTracks_;
   Float_t         TAPPionTracks_fCoordinates_fPt[kMaxTAPPionTracks];   //[TAPPionTracks_]
   Float_t         TAPPionTracks_fCoordinates_fEta[kMaxTAPPionTracks];   //[TAPPionTracks_]
   Float_t         TAPPionTracks_fCoordinates_fPhi[kMaxTAPPionTracks];   //[TAPPionTracks_]
   Float_t         TAPPionTracks_fCoordinates_fE[kMaxTAPPionTracks];   //[TAPPionTracks_]
   vector<float>   *TAPPionTracks_dxypv;
   vector<bool>    *TAPPionTracks_leptonMatch;
   vector<float>   *TAPPionTracks_mT;
   vector<float>   *TAPPionTracks_pfRelIso03chg;
   vector<float>   *TAPPionTracks_trkiso;
   vector<int>     *TriggerPass;
   vector<int>     *TriggerPrescales;
   vector<int>     *TriggerVersion;

   // List of branches
   TBranch        *b_RunNum;   //!
   TBranch        *b_LumiBlockNum;   //!
   TBranch        *b_EvtNum;   //!
   TBranch        *b_BadChargedCandidateFilter;   //!
   TBranch        *b_BadPFMuonDzFilter;   //!
   TBranch        *b_BadPFMuonFilter;   //!
   TBranch        *b_BTagsDeepCSV;   //!
   TBranch        *b_CaloMET;   //!
   TBranch        *b_CaloMETPhi;   //!
   TBranch        *b_CSCTightHaloFilter;   //!
   TBranch        *b_DeltaPhi1;   //!
   TBranch        *b_DeltaPhi1_AK8;   //!
   TBranch        *b_DeltaPhi2;   //!
   TBranch        *b_DeltaPhi2_AK8;   //!
   TBranch        *b_DeltaPhi3;   //!
   TBranch        *b_DeltaPhi4;   //!
   TBranch        *b_DeltaPhiMin_AK8;   //!
   TBranch        *b_ecalBadCalibFilter;   //!
   TBranch        *b_EcalDeadCellBoundaryEnergyFilter;   //!
   TBranch        *b_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_eeBadScFilter;   //!
   TBranch        *b_Electrons_;   //!
   TBranch        *b_Electrons_fCoordinates_fPt;   //!
   TBranch        *b_Electrons_fCoordinates_fEta;   //!
   TBranch        *b_Electrons_fCoordinates_fPhi;   //!
   TBranch        *b_Electrons_fCoordinates_fE;   //!
   TBranch        *b_Electrons_charge;   //!
   TBranch        *b_Electrons_iso;   //!
   TBranch        *b_Electrons_mediumID;   //!
   TBranch        *b_Electrons_MTW;   //!
   TBranch        *b_Electrons_passIso;   //!
   TBranch        *b_Electrons_tightID;   //!
   TBranch        *b_fixedGridRhoFastjetAll;   //!
   TBranch        *b_GenMT2_AK8;   //!
   TBranch        *b_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_globalTightHalo2016Filter;   //!
   TBranch        *b_hasGenPromptPhoton;   //!
   TBranch        *b_HBHEIsoNoiseFilter;   //!
   TBranch        *b_HBHENoiseFilter;   //!
   TBranch        *b_hfNoisyHitsFilter;   //!
   TBranch        *b_HLTMuonObjects_;   //!
   TBranch        *b_HLTMuonObjects_fCoordinates_fPt;   //!
   TBranch        *b_HLTMuonObjects_fCoordinates_fEta;   //!
   TBranch        *b_HLTMuonObjects_fCoordinates_fPhi;   //!
   TBranch        *b_HLTMuonObjects_fCoordinates_fE;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_HT5;   //!
   TBranch        *b_HTOnline;   //!
   TBranch        *b_isoElectronTracks;   //!
   TBranch        *b_isoMuonTracks;   //!
   TBranch        *b_isoPionTracks;   //!
   TBranch        *b_JetID;   //!
   TBranch        *b_JetIDAK15;   //!
   TBranch        *b_JetIDAK8;   //!
   TBranch        *b_Jets_;   //!
   TBranch        *b_Jets_fCoordinates_fPt;   //!
   TBranch        *b_Jets_fCoordinates_fEta;   //!
   TBranch        *b_Jets_fCoordinates_fPhi;   //!
   TBranch        *b_Jets_fCoordinates_fE;   //!
   TBranch        *b_Jets_axismajor;   //!
   TBranch        *b_Jets_axisminor;   //!
   TBranch        *b_Jets_bDiscriminatorCSV;   //!
   TBranch        *b_Jets_bJetTagDeepCSVBvsAll;   //!
   TBranch        *b_Jets_bJetTagDeepCSVprobb;   //!
   TBranch        *b_Jets_bJetTagDeepCSVprobbb;   //!
   TBranch        *b_Jets_bJetTagDeepCSVprobc;   //!
   TBranch        *b_Jets_bJetTagDeepCSVprobudsg;   //!
   TBranch        *b_Jets_bJetTagDeepFlavourprobb;   //!
   TBranch        *b_Jets_bJetTagDeepFlavourprobbb;   //!
   TBranch        *b_Jets_bJetTagDeepFlavourprobc;   //!
   TBranch        *b_Jets_bJetTagDeepFlavourprobg;   //!
   TBranch        *b_Jets_bJetTagDeepFlavourproblepb;   //!
   TBranch        *b_Jets_bJetTagDeepFlavourprobuds;   //!
   TBranch        *b_Jets_chargedEmEnergyFraction;   //!
   TBranch        *b_Jets_chargedHadronEnergyFraction;   //!
   TBranch        *b_Jets_chargedHadronMultiplicity;   //!
   TBranch        *b_Jets_chargedMultiplicity;   //!
   TBranch        *b_Jets_electronEnergyFraction;   //!
   TBranch        *b_Jets_electronMultiplicity;   //!
   TBranch        *b_Jets_hadronFlavor;   //!
   TBranch        *b_Jets_hfEMEnergyFraction;   //!
   TBranch        *b_Jets_hfHadronEnergyFraction;   //!
   TBranch        *b_Jets_HTMask;   //!
   TBranch        *b_Jets_ID;   //!
   TBranch        *b_Jets_jecFactor;   //!
   TBranch        *b_Jets_jecUnc;   //!
   TBranch        *b_Jets_LeptonMask;   //!
   TBranch        *b_Jets_MHTMask;   //!
   TBranch        *b_Jets_multiplicity;   //!
   TBranch        *b_Jets_muonEnergyFraction;   //!
   TBranch        *b_Jets_muonMultiplicity;   //!
   TBranch        *b_Jets_neutralEmEnergyFraction;   //!
   TBranch        *b_Jets_neutralHadronEnergyFraction;   //!
   TBranch        *b_Jets_neutralHadronMultiplicity;   //!
   TBranch        *b_Jets_neutralMultiplicity;   //!
   TBranch        *b_Jets_partonFlavor;   //!
   TBranch        *b_Jets_photonEnergyFraction;   //!
   TBranch        *b_Jets_photonMultiplicity;   //!
   TBranch        *b_Jets_pileupId;   //!
   TBranch        *b_Jets_ptD;   //!
   TBranch        *b_Jets_qgLikelihood;   //!
   TBranch        *b_JetsAK15_;   //!
   TBranch        *b_JetsAK15_fCoordinates_fPt;   //!
   TBranch        *b_JetsAK15_fCoordinates_fEta;   //!
   TBranch        *b_JetsAK15_fCoordinates_fPhi;   //!
   TBranch        *b_JetsAK15_fCoordinates_fE;   //!
   TBranch        *b_JetsAK15_axismajor;   //!
   TBranch        *b_JetsAK15_axisminor;   //!
   TBranch        *b_JetsAK15_chargedEmEnergyFraction;   //!
   TBranch        *b_JetsAK15_chargedHadronEnergyFraction;   //!
   TBranch        *b_JetsAK15_chargedHadronMultiplicity;   //!
   TBranch        *b_JetsAK15_chargedMultiplicity;   //!
   TBranch        *b_JetsAK15_constituentsIndex;   //!
   TBranch        *b_JetsAK15_constituentsIndexCounts;   //!
   TBranch        *b_JetsAK15_DeepMassDecorrelTagbbvsLight;   //!
   TBranch        *b_JetsAK15_DeepMassDecorrelTagHbbvsQCD;   //!
   TBranch        *b_JetsAK15_DeepMassDecorrelTagTvsQCD;   //!
   TBranch        *b_JetsAK15_DeepMassDecorrelTagWvsQCD;   //!
   TBranch        *b_JetsAK15_DeepMassDecorrelTagZbbvsQCD;   //!
   TBranch        *b_JetsAK15_DeepMassDecorrelTagZHbbvsQCD;   //!
   TBranch        *b_JetsAK15_DeepMassDecorrelTagZvsQCD;   //!
   TBranch        *b_JetsAK15_DeepTagHbbvsQCD;   //!
   TBranch        *b_JetsAK15_DeepTagTvsQCD;   //!
   TBranch        *b_JetsAK15_DeepTagWvsQCD;   //!
   TBranch        *b_JetsAK15_DeepTagZbbvsQCD;   //!
   TBranch        *b_JetsAK15_DeepTagZvsQCD;   //!
   TBranch        *b_JetsAK15_doubleBDiscriminator;   //!
   TBranch        *b_JetsAK15_ecfC2b1;   //!
   TBranch        *b_JetsAK15_ecfC2b2;   //!
   TBranch        *b_JetsAK15_ecfD2b1;   //!
   TBranch        *b_JetsAK15_ecfD2b2;   //!
   TBranch        *b_JetsAK15_ecfM2b1;   //!
   TBranch        *b_JetsAK15_ecfM2b2;   //!
   TBranch        *b_JetsAK15_ecfN2b1;   //!
   TBranch        *b_JetsAK15_ecfN2b2;   //!
   TBranch        *b_JetsAK15_electronEnergyFraction;   //!
   TBranch        *b_JetsAK15_electronMultiplicity;   //!
   TBranch        *b_JetsAK15_girth;   //!
   TBranch        *b_JetsAK15_hadronFlavor;   //!
   TBranch        *b_JetsAK15_hfEMEnergyFraction;   //!
   TBranch        *b_JetsAK15_hfHadronEnergyFraction;   //!
   TBranch        *b_JetsAK15_ID;   //!
   TBranch        *b_JetsAK15_jecFactor;   //!
   TBranch        *b_JetsAK15_multiplicity;   //!
   TBranch        *b_JetsAK15_muonEnergyFraction;   //!
   TBranch        *b_JetsAK15_muonMultiplicity;   //!
   TBranch        *b_JetsAK15_neutralEmEnergyFraction;   //!
   TBranch        *b_JetsAK15_neutralHadronEnergyFraction;   //!
   TBranch        *b_JetsAK15_neutralHadronMultiplicity;   //!
   TBranch        *b_JetsAK15_neutralMultiplicity;   //!
   TBranch        *b_JetsAK15_NsubjettinessTau1;   //!
   TBranch        *b_JetsAK15_NsubjettinessTau2;   //!
   TBranch        *b_JetsAK15_NsubjettinessTau3;   //!
   TBranch        *b_JetsAK15_NsubjettinessTau4;   //!
   TBranch        *b_JetsAK15_NumBhadrons;   //!
   TBranch        *b_JetsAK15_NumChadrons;   //!
   TBranch        *b_JetsAK15_partonFlavor;   //!
   TBranch        *b_JetsAK15_pfMassIndependentDeepDoubleBvLJetTagsProbHbb;   //!
   TBranch        *b_JetsAK15_photonEnergyFraction;   //!
   TBranch        *b_JetsAK15_photonMultiplicity;   //!
   TBranch        *b_JetsAK15_ptD;   //!
   TBranch        *b_JetsAK15_softDropMass;   //!
   TBranch        *b_JetsAK15_softDropMassBeta1;   //!
   TBranch        *b_JetsAK15_subjets_;   //!
   TBranch        *b_JetsAK15_subjets_fCoordinates_fPt;   //!
   TBranch        *b_JetsAK15_subjets_fCoordinates_fEta;   //!
   TBranch        *b_JetsAK15_subjets_fCoordinates_fPhi;   //!
   TBranch        *b_JetsAK15_subjets_fCoordinates_fE;   //!
   TBranch        *b_JetsAK15_subjetsCounts;   //!
   TBranch        *b_JetsAK8_;   //!
   TBranch        *b_JetsAK8_fCoordinates_fPt;   //!
   TBranch        *b_JetsAK8_fCoordinates_fEta;   //!
   TBranch        *b_JetsAK8_fCoordinates_fPhi;   //!
   TBranch        *b_JetsAK8_fCoordinates_fE;   //!
   TBranch        *b_JetsAK8_axismajor;   //!
   TBranch        *b_JetsAK8_axisminor;   //!
   TBranch        *b_JetsAK8_chargedEmEnergyFraction;   //!
   TBranch        *b_JetsAK8_chargedHadronEnergyFraction;   //!
   TBranch        *b_JetsAK8_chargedHadronMultiplicity;   //!
   TBranch        *b_JetsAK8_chargedMultiplicity;   //!
   TBranch        *b_JetsAK8_constituentsIndex;   //!
   TBranch        *b_JetsAK8_constituentsIndexCounts;   //!
   TBranch        *b_JetsAK8_DeepMassDecorrelTagbbvsLight;   //!
   TBranch        *b_JetsAK8_DeepMassDecorrelTagHbbvsQCD;   //!
   TBranch        *b_JetsAK8_DeepMassDecorrelTagTvsQCD;   //!
   TBranch        *b_JetsAK8_DeepMassDecorrelTagWvsQCD;   //!
   TBranch        *b_JetsAK8_DeepMassDecorrelTagZbbvsQCD;   //!
   TBranch        *b_JetsAK8_DeepMassDecorrelTagZHbbvsQCD;   //!
   TBranch        *b_JetsAK8_DeepMassDecorrelTagZvsQCD;   //!
   TBranch        *b_JetsAK8_DeepTagHbbvsQCD;   //!
   TBranch        *b_JetsAK8_DeepTagTvsQCD;   //!
   TBranch        *b_JetsAK8_DeepTagWvsQCD;   //!
   TBranch        *b_JetsAK8_DeepTagZbbvsQCD;   //!
   TBranch        *b_JetsAK8_DeepTagZvsQCD;   //!
   TBranch        *b_JetsAK8_doubleBDiscriminator;   //!
   TBranch        *b_JetsAK8_ecfN2b1;   //!
   TBranch        *b_JetsAK8_ecfN2b2;   //!
   TBranch        *b_JetsAK8_ecfN3b1;   //!
   TBranch        *b_JetsAK8_ecfN3b2;   //!
   TBranch        *b_JetsAK8_electronEnergyFraction;   //!
   TBranch        *b_JetsAK8_electronMultiplicity;   //!
   TBranch        *b_JetsAK8_girth;   //!
   TBranch        *b_JetsAK8_hadronFlavor;   //!
   TBranch        *b_JetsAK8_hfEMEnergyFraction;   //!
   TBranch        *b_JetsAK8_hfHadronEnergyFraction;   //!
   TBranch        *b_JetsAK8_ID;   //!
   TBranch        *b_JetsAK8_jecFactor;   //!
   TBranch        *b_JetsAK8_jecUnc;   //!
   TBranch        *b_JetsAK8_multiplicity;   //!
   TBranch        *b_JetsAK8_muonEnergyFraction;   //!
   TBranch        *b_JetsAK8_muonMultiplicity;   //!
   TBranch        *b_JetsAK8_neutralEmEnergyFraction;   //!
   TBranch        *b_JetsAK8_neutralHadronEnergyFraction;   //!
   TBranch        *b_JetsAK8_neutralHadronMultiplicity;   //!
   TBranch        *b_JetsAK8_neutralMultiplicity;   //!
   TBranch        *b_JetsAK8_NsubjettinessTau1;   //!
   TBranch        *b_JetsAK8_NsubjettinessTau2;   //!
   TBranch        *b_JetsAK8_NsubjettinessTau3;   //!
   TBranch        *b_JetsAK8_NumBhadrons;   //!
   TBranch        *b_JetsAK8_NumChadrons;   //!
   TBranch        *b_JetsAK8_partonFlavor;   //!
   TBranch        *b_JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb;   //!
   TBranch        *b_JetsAK8_photonEnergyFraction;   //!
   TBranch        *b_JetsAK8_photonMultiplicity;   //!
   TBranch        *b_JetsAK8_ptD;   //!
   TBranch        *b_JetsAK8_softDropMass;   //!
   TBranch        *b_JetsAK8_subjets_;   //!
   TBranch        *b_JetsAK8_subjets_fCoordinates_fPt;   //!
   TBranch        *b_JetsAK8_subjets_fCoordinates_fEta;   //!
   TBranch        *b_JetsAK8_subjets_fCoordinates_fPhi;   //!
   TBranch        *b_JetsAK8_subjets_fCoordinates_fE;   //!
   TBranch        *b_JetsAK8_subjetsCounts;   //!
   TBranch        *b_JetsAK8_subjets_axismajor;   //!
   TBranch        *b_JetsAK8_subjets_axisminor;   //!
   TBranch        *b_JetsAK8_subjets_jecFactor;   //!
   TBranch        *b_JetsAK8_subjets_multiplicity;   //!
   TBranch        *b_JetsAK8_subjets_ptD;   //!
   TBranch        *b_JetsConstituents_;   //!
   TBranch        *b_JetsConstituents_fCoordinates_fPt;   //!
   TBranch        *b_JetsConstituents_fCoordinates_fEta;   //!
   TBranch        *b_JetsConstituents_fCoordinates_fPhi;   //!
   TBranch        *b_JetsConstituents_fCoordinates_fE;   //!
   TBranch        *b_JetsConstituents_dxy;   //!
   TBranch        *b_JetsConstituents_dxysig;   //!
   TBranch        *b_JetsConstituents_dz;   //!
   TBranch        *b_JetsConstituents_dzsig;   //!
   TBranch        *b_JetsConstituents_PdgId;   //!
   TBranch        *b_JetsConstituents_PuppiWeight;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METPhi;   //!
   TBranch        *b_METSignificance;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_MHTOnline;   //!
   TBranch        *b_MHTPhi;   //!
   TBranch        *b_MJJ_AK8;   //!
   TBranch        *b_Mmc_AK8;   //!
   TBranch        *b_MT_AK8;   //!
   TBranch        *b_Muons_;   //!
   TBranch        *b_Muons_fCoordinates_fPt;   //!
   TBranch        *b_Muons_fCoordinates_fEta;   //!
   TBranch        *b_Muons_fCoordinates_fPhi;   //!
   TBranch        *b_Muons_fCoordinates_fE;   //!
   TBranch        *b_Muons_charge;   //!
   TBranch        *b_Muons_iso;   //!
   TBranch        *b_Muons_mediumID;   //!
   TBranch        *b_Muons_MTW;   //!
   TBranch        *b_Muons_passIso;   //!
   TBranch        *b_Muons_tightID;   //!
   TBranch        *b_nAllVertices;   //!
   TBranch        *b_NElectrons;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_NMuons;   //!
   TBranch        *b_NonPrefiringProb;   //!
   TBranch        *b_NonPrefiringProbDown;   //!
   TBranch        *b_NonPrefiringProbECAL;   //!
   TBranch        *b_NonPrefiringProbECALDown;   //!
   TBranch        *b_NonPrefiringProbECALUp;   //!
   TBranch        *b_NonPrefiringProbMuon;   //!
   TBranch        *b_NonPrefiringProbMuonDown;   //!
   TBranch        *b_NonPrefiringProbMuonUp;   //!
   TBranch        *b_NonPrefiringProbUp;   //!
   TBranch        *b_NVtx;   //!
   TBranch        *b_PFCaloMETRatio;   //!
   TBranch        *b_Photons_;   //!
   TBranch        *b_Photons_fCoordinates_fPt;   //!
   TBranch        *b_Photons_fCoordinates_fEta;   //!
   TBranch        *b_Photons_fCoordinates_fPhi;   //!
   TBranch        *b_Photons_fCoordinates_fE;   //!
   TBranch        *b_Photons_electronFakes;   //!
   TBranch        *b_Photons_fullID;   //!
   TBranch        *b_Photons_genMatched;   //!
   TBranch        *b_Photons_hadTowOverEM;   //!
   TBranch        *b_Photons_hasPixelSeed;   //!
   TBranch        *b_Photons_isEB;   //!
   TBranch        *b_Photons_nonPrompt;   //!
   TBranch        *b_Photons_passElectronVeto;   //!
   TBranch        *b_Photons_pfChargedIso;   //!
   TBranch        *b_Photons_pfChargedIsoRhoCorr;   //!
   TBranch        *b_Photons_pfGammaIso;   //!
   TBranch        *b_Photons_pfGammaIsoRhoCorr;   //!
   TBranch        *b_Photons_pfNeutralIso;   //!
   TBranch        *b_Photons_pfNeutralIsoRhoCorr;   //!
   TBranch        *b_Photons_sigmaIetaIeta;   //!
   TBranch        *b_PrescaleWeightHT;   //!
   TBranch        *b_PrimaryVertexFilter;   //!
   TBranch        *b_TAPElectronTracks_;   //!
   TBranch        *b_TAPElectronTracks_fCoordinates_fPt;   //!
   TBranch        *b_TAPElectronTracks_fCoordinates_fEta;   //!
   TBranch        *b_TAPElectronTracks_fCoordinates_fPhi;   //!
   TBranch        *b_TAPElectronTracks_fCoordinates_fE;   //!
   TBranch        *b_TAPElectronTracks_dxypv;   //!
   TBranch        *b_TAPElectronTracks_leptonMatch;   //!
   TBranch        *b_TAPElectronTracks_mT;   //!
   TBranch        *b_TAPElectronTracks_pfRelIso03chg;   //!
   TBranch        *b_TAPElectronTracks_trkiso;   //!
   TBranch        *b_TAPMuonTracks_;   //!
   TBranch        *b_TAPMuonTracks_fCoordinates_fPt;   //!
   TBranch        *b_TAPMuonTracks_fCoordinates_fEta;   //!
   TBranch        *b_TAPMuonTracks_fCoordinates_fPhi;   //!
   TBranch        *b_TAPMuonTracks_fCoordinates_fE;   //!
   TBranch        *b_TAPMuonTracks_dxypv;   //!
   TBranch        *b_TAPMuonTracks_leptonMatch;   //!
   TBranch        *b_TAPMuonTracks_mT;   //!
   TBranch        *b_TAPMuonTracks_pfRelIso03chg;   //!
   TBranch        *b_TAPMuonTracks_trkiso;   //!
   TBranch        *b_TAPPionTracks_;   //!
   TBranch        *b_TAPPionTracks_fCoordinates_fPt;   //!
   TBranch        *b_TAPPionTracks_fCoordinates_fEta;   //!
   TBranch        *b_TAPPionTracks_fCoordinates_fPhi;   //!
   TBranch        *b_TAPPionTracks_fCoordinates_fE;   //!
   TBranch        *b_TAPPionTracks_dxypv;   //!
   TBranch        *b_TAPPionTracks_leptonMatch;   //!
   TBranch        *b_TAPPionTracks_mT;   //!
   TBranch        *b_TAPPionTracks_pfRelIso03chg;   //!
   TBranch        *b_TAPPionTracks_trkiso;   //!
   TBranch        *b_TriggerPass;   //!
   TBranch        *b_TriggerPrescales;   //!
   TBranch        *b_TriggerVersion;   //!

   allin1(TTree *tree=0);
   virtual ~allin1();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef allin1_cxx
allin1::allin1(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2018A-UL2018-v3/SingleMuon/2585_RA2AnalysisTree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2018A-UL2018-v3/SingleMuon/2585_RA2AnalysisTree.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2018A-UL2018-v3/SingleMuon/2585_RA2AnalysisTree.root:/TreeMaker2");
      dir->GetObject("PreSelection",tree);

   }
   Init(tree);
}

allin1::~allin1()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t allin1::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t allin1::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void allin1::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Electrons_charge = 0;
   Electrons_iso = 0;
   Electrons_mediumID = 0;
   Electrons_MTW = 0;
   Electrons_passIso = 0;
   Electrons_tightID = 0;
   Jets_axismajor = 0;
   Jets_axisminor = 0;
   Jets_bDiscriminatorCSV = 0;
   Jets_bJetTagDeepCSVBvsAll = 0;
   Jets_bJetTagDeepCSVprobb = 0;
   Jets_bJetTagDeepCSVprobbb = 0;
   Jets_bJetTagDeepCSVprobc = 0;
   Jets_bJetTagDeepCSVprobudsg = 0;
   Jets_bJetTagDeepFlavourprobb = 0;
   Jets_bJetTagDeepFlavourprobbb = 0;
   Jets_bJetTagDeepFlavourprobc = 0;
   Jets_bJetTagDeepFlavourprobg = 0;
   Jets_bJetTagDeepFlavourproblepb = 0;
   Jets_bJetTagDeepFlavourprobuds = 0;
   Jets_chargedEmEnergyFraction = 0;
   Jets_chargedHadronEnergyFraction = 0;
   Jets_chargedHadronMultiplicity = 0;
   Jets_chargedMultiplicity = 0;
   Jets_electronEnergyFraction = 0;
   Jets_electronMultiplicity = 0;
   Jets_hadronFlavor = 0;
   Jets_hfEMEnergyFraction = 0;
   Jets_hfHadronEnergyFraction = 0;
   Jets_HTMask = 0;
   Jets_ID = 0;
   Jets_jecFactor = 0;
   Jets_jecUnc = 0;
   Jets_LeptonMask = 0;
   Jets_MHTMask = 0;
   Jets_multiplicity = 0;
   Jets_muonEnergyFraction = 0;
   Jets_muonMultiplicity = 0;
   Jets_neutralEmEnergyFraction = 0;
   Jets_neutralHadronEnergyFraction = 0;
   Jets_neutralHadronMultiplicity = 0;
   Jets_neutralMultiplicity = 0;
   Jets_partonFlavor = 0;
   Jets_photonEnergyFraction = 0;
   Jets_photonMultiplicity = 0;
   Jets_pileupId = 0;
   Jets_ptD = 0;
   Jets_qgLikelihood = 0;
   JetsAK15_axismajor = 0;
   JetsAK15_axisminor = 0;
   JetsAK15_chargedEmEnergyFraction = 0;
   JetsAK15_chargedHadronEnergyFraction = 0;
   JetsAK15_chargedHadronMultiplicity = 0;
   JetsAK15_chargedMultiplicity = 0;
   JetsAK15_constituentsIndex = 0;
   JetsAK15_constituentsIndexCounts = 0;
   JetsAK15_DeepMassDecorrelTagbbvsLight = 0;
   JetsAK15_DeepMassDecorrelTagHbbvsQCD = 0;
   JetsAK15_DeepMassDecorrelTagTvsQCD = 0;
   JetsAK15_DeepMassDecorrelTagWvsQCD = 0;
   JetsAK15_DeepMassDecorrelTagZbbvsQCD = 0;
   JetsAK15_DeepMassDecorrelTagZHbbvsQCD = 0;
   JetsAK15_DeepMassDecorrelTagZvsQCD = 0;
   JetsAK15_DeepTagHbbvsQCD = 0;
   JetsAK15_DeepTagTvsQCD = 0;
   JetsAK15_DeepTagWvsQCD = 0;
   JetsAK15_DeepTagZbbvsQCD = 0;
   JetsAK15_DeepTagZvsQCD = 0;
   JetsAK15_doubleBDiscriminator = 0;
   JetsAK15_ecfC2b1 = 0;
   JetsAK15_ecfC2b2 = 0;
   JetsAK15_ecfD2b1 = 0;
   JetsAK15_ecfD2b2 = 0;
   JetsAK15_ecfM2b1 = 0;
   JetsAK15_ecfM2b2 = 0;
   JetsAK15_ecfN2b1 = 0;
   JetsAK15_ecfN2b2 = 0;
   JetsAK15_electronEnergyFraction = 0;
   JetsAK15_electronMultiplicity = 0;
   JetsAK15_girth = 0;
   JetsAK15_hadronFlavor = 0;
   JetsAK15_hfEMEnergyFraction = 0;
   JetsAK15_hfHadronEnergyFraction = 0;
   JetsAK15_ID = 0;
   JetsAK15_jecFactor = 0;
   JetsAK15_multiplicity = 0;
   JetsAK15_muonEnergyFraction = 0;
   JetsAK15_muonMultiplicity = 0;
   JetsAK15_neutralEmEnergyFraction = 0;
   JetsAK15_neutralHadronEnergyFraction = 0;
   JetsAK15_neutralHadronMultiplicity = 0;
   JetsAK15_neutralMultiplicity = 0;
   JetsAK15_NsubjettinessTau1 = 0;
   JetsAK15_NsubjettinessTau2 = 0;
   JetsAK15_NsubjettinessTau3 = 0;
   JetsAK15_NsubjettinessTau4 = 0;
   JetsAK15_NumBhadrons = 0;
   JetsAK15_NumChadrons = 0;
   JetsAK15_partonFlavor = 0;
   JetsAK15_pfMassIndependentDeepDoubleBvLJetTagsProbHbb = 0;
   JetsAK15_photonEnergyFraction = 0;
   JetsAK15_photonMultiplicity = 0;
   JetsAK15_ptD = 0;
   JetsAK15_softDropMass = 0;
   JetsAK15_softDropMassBeta1 = 0;
   JetsAK15_subjetsCounts = 0;
   JetsAK8_axismajor = 0;
   JetsAK8_axisminor = 0;
   JetsAK8_chargedEmEnergyFraction = 0;
   JetsAK8_chargedHadronEnergyFraction = 0;
   JetsAK8_chargedHadronMultiplicity = 0;
   JetsAK8_chargedMultiplicity = 0;
   JetsAK8_constituentsIndex = 0;
   JetsAK8_constituentsIndexCounts = 0;
   JetsAK8_DeepMassDecorrelTagbbvsLight = 0;
   JetsAK8_DeepMassDecorrelTagHbbvsQCD = 0;
   JetsAK8_DeepMassDecorrelTagTvsQCD = 0;
   JetsAK8_DeepMassDecorrelTagWvsQCD = 0;
   JetsAK8_DeepMassDecorrelTagZbbvsQCD = 0;
   JetsAK8_DeepMassDecorrelTagZHbbvsQCD = 0;
   JetsAK8_DeepMassDecorrelTagZvsQCD = 0;
   JetsAK8_DeepTagHbbvsQCD = 0;
   JetsAK8_DeepTagTvsQCD = 0;
   JetsAK8_DeepTagWvsQCD = 0;
   JetsAK8_DeepTagZbbvsQCD = 0;
   JetsAK8_DeepTagZvsQCD = 0;
   JetsAK8_doubleBDiscriminator = 0;
   JetsAK8_ecfN2b1 = 0;
   JetsAK8_ecfN2b2 = 0;
   JetsAK8_ecfN3b1 = 0;
   JetsAK8_ecfN3b2 = 0;
   JetsAK8_electronEnergyFraction = 0;
   JetsAK8_electronMultiplicity = 0;
   JetsAK8_girth = 0;
   JetsAK8_hadronFlavor = 0;
   JetsAK8_hfEMEnergyFraction = 0;
   JetsAK8_hfHadronEnergyFraction = 0;
   JetsAK8_ID = 0;
   JetsAK8_jecFactor = 0;
   JetsAK8_jecUnc = 0;
   JetsAK8_multiplicity = 0;
   JetsAK8_muonEnergyFraction = 0;
   JetsAK8_muonMultiplicity = 0;
   JetsAK8_neutralEmEnergyFraction = 0;
   JetsAK8_neutralHadronEnergyFraction = 0;
   JetsAK8_neutralHadronMultiplicity = 0;
   JetsAK8_neutralMultiplicity = 0;
   JetsAK8_NsubjettinessTau1 = 0;
   JetsAK8_NsubjettinessTau2 = 0;
   JetsAK8_NsubjettinessTau3 = 0;
   JetsAK8_NumBhadrons = 0;
   JetsAK8_NumChadrons = 0;
   JetsAK8_partonFlavor = 0;
   JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb = 0;
   JetsAK8_photonEnergyFraction = 0;
   JetsAK8_photonMultiplicity = 0;
   JetsAK8_ptD = 0;
   JetsAK8_softDropMass = 0;
   JetsAK8_subjetsCounts = 0;
   JetsAK8_subjets_axismajor = 0;
   JetsAK8_subjets_axisminor = 0;
   JetsAK8_subjets_jecFactor = 0;
   JetsAK8_subjets_multiplicity = 0;
   JetsAK8_subjets_ptD = 0;
   JetsConstituents_dxy = 0;
   JetsConstituents_dxysig = 0;
   JetsConstituents_dz = 0;
   JetsConstituents_dzsig = 0;
   JetsConstituents_PdgId = 0;
   JetsConstituents_PuppiWeight = 0;
   Muons_charge = 0;
   Muons_iso = 0;
   Muons_mediumID = 0;
   Muons_MTW = 0;
   Muons_passIso = 0;
   Muons_tightID = 0;
   Photons_electronFakes = 0;
   Photons_fullID = 0;
   Photons_genMatched = 0;
   Photons_hadTowOverEM = 0;
   Photons_hasPixelSeed = 0;
   Photons_isEB = 0;
   Photons_nonPrompt = 0;
   Photons_passElectronVeto = 0;
   Photons_pfChargedIso = 0;
   Photons_pfChargedIsoRhoCorr = 0;
   Photons_pfGammaIso = 0;
   Photons_pfGammaIsoRhoCorr = 0;
   Photons_pfNeutralIso = 0;
   Photons_pfNeutralIsoRhoCorr = 0;
   Photons_sigmaIetaIeta = 0;
   TAPElectronTracks_dxypv = 0;
   TAPElectronTracks_leptonMatch = 0;
   TAPElectronTracks_mT = 0;
   TAPElectronTracks_pfRelIso03chg = 0;
   TAPElectronTracks_trkiso = 0;
   TAPMuonTracks_dxypv = 0;
   TAPMuonTracks_leptonMatch = 0;
   TAPMuonTracks_mT = 0;
   TAPMuonTracks_pfRelIso03chg = 0;
   TAPMuonTracks_trkiso = 0;
   TAPPionTracks_dxypv = 0;
   TAPPionTracks_leptonMatch = 0;
   TAPPionTracks_mT = 0;
   TAPPionTracks_pfRelIso03chg = 0;
   TAPPionTracks_trkiso = 0;
   TriggerPass = 0;
   TriggerPrescales = 0;
   TriggerVersion = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNum", &RunNum, &b_RunNum);
   fChain->SetBranchAddress("LumiBlockNum", &LumiBlockNum, &b_LumiBlockNum);
   fChain->SetBranchAddress("EvtNum", &EvtNum, &b_EvtNum);
   fChain->SetBranchAddress("BadChargedCandidateFilter", &BadChargedCandidateFilter, &b_BadChargedCandidateFilter);
   fChain->SetBranchAddress("BadPFMuonDzFilter", &BadPFMuonDzFilter, &b_BadPFMuonDzFilter);
   fChain->SetBranchAddress("BadPFMuonFilter", &BadPFMuonFilter, &b_BadPFMuonFilter);
   fChain->SetBranchAddress("BTagsDeepCSV", &BTagsDeepCSV, &b_BTagsDeepCSV);
   fChain->SetBranchAddress("CaloMET", &CaloMET, &b_CaloMET);
   fChain->SetBranchAddress("CaloMETPhi", &CaloMETPhi, &b_CaloMETPhi);
   fChain->SetBranchAddress("CSCTightHaloFilter", &CSCTightHaloFilter, &b_CSCTightHaloFilter);
   fChain->SetBranchAddress("DeltaPhi1", &DeltaPhi1, &b_DeltaPhi1);
   fChain->SetBranchAddress("DeltaPhi1_AK8", &DeltaPhi1_AK8, &b_DeltaPhi1_AK8);
   fChain->SetBranchAddress("DeltaPhi2", &DeltaPhi2, &b_DeltaPhi2);
   fChain->SetBranchAddress("DeltaPhi2_AK8", &DeltaPhi2_AK8, &b_DeltaPhi2_AK8);
   fChain->SetBranchAddress("DeltaPhi3", &DeltaPhi3, &b_DeltaPhi3);
   fChain->SetBranchAddress("DeltaPhi4", &DeltaPhi4, &b_DeltaPhi4);
   fChain->SetBranchAddress("DeltaPhiMin_AK8", &DeltaPhiMin_AK8, &b_DeltaPhiMin_AK8);
   fChain->SetBranchAddress("ecalBadCalibFilter", &ecalBadCalibFilter, &b_ecalBadCalibFilter);
   fChain->SetBranchAddress("EcalDeadCellBoundaryEnergyFilter", &EcalDeadCellBoundaryEnergyFilter, &b_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("EcalDeadCellTriggerPrimitiveFilter", &EcalDeadCellTriggerPrimitiveFilter, &b_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("eeBadScFilter", &eeBadScFilter, &b_eeBadScFilter);
   fChain->SetBranchAddress("Electrons", &Electrons_, &b_Electrons_);
   fChain->SetBranchAddress("Electrons.fCoordinates.fPt", Electrons_fCoordinates_fPt, &b_Electrons_fCoordinates_fPt);
   fChain->SetBranchAddress("Electrons.fCoordinates.fEta", Electrons_fCoordinates_fEta, &b_Electrons_fCoordinates_fEta);
   fChain->SetBranchAddress("Electrons.fCoordinates.fPhi", Electrons_fCoordinates_fPhi, &b_Electrons_fCoordinates_fPhi);
   fChain->SetBranchAddress("Electrons.fCoordinates.fE", Electrons_fCoordinates_fE, &b_Electrons_fCoordinates_fE);
   fChain->SetBranchAddress("Electrons_charge", &Electrons_charge, &b_Electrons_charge);
   fChain->SetBranchAddress("Electrons_iso", &Electrons_iso, &b_Electrons_iso);
   fChain->SetBranchAddress("Electrons_mediumID", &Electrons_mediumID, &b_Electrons_mediumID);
   fChain->SetBranchAddress("Electrons_MTW", &Electrons_MTW, &b_Electrons_MTW);
   fChain->SetBranchAddress("Electrons_passIso", &Electrons_passIso, &b_Electrons_passIso);
   fChain->SetBranchAddress("Electrons_tightID", &Electrons_tightID, &b_Electrons_tightID);
   fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("GenMT2_AK8", &GenMT2_AK8, &b_GenMT2_AK8);
   fChain->SetBranchAddress("globalSuperTightHalo2016Filter", &globalSuperTightHalo2016Filter, &b_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("globalTightHalo2016Filter", &globalTightHalo2016Filter, &b_globalTightHalo2016Filter);
   fChain->SetBranchAddress("hasGenPromptPhoton", &hasGenPromptPhoton, &b_hasGenPromptPhoton);
   fChain->SetBranchAddress("HBHEIsoNoiseFilter", &HBHEIsoNoiseFilter, &b_HBHEIsoNoiseFilter);
   fChain->SetBranchAddress("HBHENoiseFilter", &HBHENoiseFilter, &b_HBHENoiseFilter);
   fChain->SetBranchAddress("hfNoisyHitsFilter", &hfNoisyHitsFilter, &b_hfNoisyHitsFilter);
   fChain->SetBranchAddress("HLTMuonObjects", &HLTMuonObjects_, &b_HLTMuonObjects_);
   fChain->SetBranchAddress("HLTMuonObjects.fCoordinates.fPt", HLTMuonObjects_fCoordinates_fPt, &b_HLTMuonObjects_fCoordinates_fPt);
   fChain->SetBranchAddress("HLTMuonObjects.fCoordinates.fEta", HLTMuonObjects_fCoordinates_fEta, &b_HLTMuonObjects_fCoordinates_fEta);
   fChain->SetBranchAddress("HLTMuonObjects.fCoordinates.fPhi", HLTMuonObjects_fCoordinates_fPhi, &b_HLTMuonObjects_fCoordinates_fPhi);
   fChain->SetBranchAddress("HLTMuonObjects.fCoordinates.fE", HLTMuonObjects_fCoordinates_fE, &b_HLTMuonObjects_fCoordinates_fE);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("HT5", &HT5, &b_HT5);
   fChain->SetBranchAddress("HTOnline", &HTOnline, &b_HTOnline);
   fChain->SetBranchAddress("isoElectronTracks", &isoElectronTracks, &b_isoElectronTracks);
   fChain->SetBranchAddress("isoMuonTracks", &isoMuonTracks, &b_isoMuonTracks);
   fChain->SetBranchAddress("isoPionTracks", &isoPionTracks, &b_isoPionTracks);
   fChain->SetBranchAddress("JetID", &JetID, &b_JetID);
   fChain->SetBranchAddress("JetIDAK15", &JetIDAK15, &b_JetIDAK15);
   fChain->SetBranchAddress("JetIDAK8", &JetIDAK8, &b_JetIDAK8);
   fChain->SetBranchAddress("Jets", &Jets_, &b_Jets_);
   fChain->SetBranchAddress("Jets.fCoordinates.fPt", Jets_fCoordinates_fPt, &b_Jets_fCoordinates_fPt);
   fChain->SetBranchAddress("Jets.fCoordinates.fEta", Jets_fCoordinates_fEta, &b_Jets_fCoordinates_fEta);
   fChain->SetBranchAddress("Jets.fCoordinates.fPhi", Jets_fCoordinates_fPhi, &b_Jets_fCoordinates_fPhi);
   fChain->SetBranchAddress("Jets.fCoordinates.fE", Jets_fCoordinates_fE, &b_Jets_fCoordinates_fE);
   fChain->SetBranchAddress("Jets_axismajor", &Jets_axismajor, &b_Jets_axismajor);
   fChain->SetBranchAddress("Jets_axisminor", &Jets_axisminor, &b_Jets_axisminor);
   fChain->SetBranchAddress("Jets_bDiscriminatorCSV", &Jets_bDiscriminatorCSV, &b_Jets_bDiscriminatorCSV);
   fChain->SetBranchAddress("Jets_bJetTagDeepCSVBvsAll", &Jets_bJetTagDeepCSVBvsAll, &b_Jets_bJetTagDeepCSVBvsAll);
   fChain->SetBranchAddress("Jets_bJetTagDeepCSVprobb", &Jets_bJetTagDeepCSVprobb, &b_Jets_bJetTagDeepCSVprobb);
   fChain->SetBranchAddress("Jets_bJetTagDeepCSVprobbb", &Jets_bJetTagDeepCSVprobbb, &b_Jets_bJetTagDeepCSVprobbb);
   fChain->SetBranchAddress("Jets_bJetTagDeepCSVprobc", &Jets_bJetTagDeepCSVprobc, &b_Jets_bJetTagDeepCSVprobc);
   fChain->SetBranchAddress("Jets_bJetTagDeepCSVprobudsg", &Jets_bJetTagDeepCSVprobudsg, &b_Jets_bJetTagDeepCSVprobudsg);
   fChain->SetBranchAddress("Jets_bJetTagDeepFlavourprobb", &Jets_bJetTagDeepFlavourprobb, &b_Jets_bJetTagDeepFlavourprobb);
   fChain->SetBranchAddress("Jets_bJetTagDeepFlavourprobbb", &Jets_bJetTagDeepFlavourprobbb, &b_Jets_bJetTagDeepFlavourprobbb);
   fChain->SetBranchAddress("Jets_bJetTagDeepFlavourprobc", &Jets_bJetTagDeepFlavourprobc, &b_Jets_bJetTagDeepFlavourprobc);
   fChain->SetBranchAddress("Jets_bJetTagDeepFlavourprobg", &Jets_bJetTagDeepFlavourprobg, &b_Jets_bJetTagDeepFlavourprobg);
   fChain->SetBranchAddress("Jets_bJetTagDeepFlavourproblepb", &Jets_bJetTagDeepFlavourproblepb, &b_Jets_bJetTagDeepFlavourproblepb);
   fChain->SetBranchAddress("Jets_bJetTagDeepFlavourprobuds", &Jets_bJetTagDeepFlavourprobuds, &b_Jets_bJetTagDeepFlavourprobuds);
   fChain->SetBranchAddress("Jets_chargedEmEnergyFraction", &Jets_chargedEmEnergyFraction, &b_Jets_chargedEmEnergyFraction);
   fChain->SetBranchAddress("Jets_chargedHadronEnergyFraction", &Jets_chargedHadronEnergyFraction, &b_Jets_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("Jets_chargedHadronMultiplicity", &Jets_chargedHadronMultiplicity, &b_Jets_chargedHadronMultiplicity);
   fChain->SetBranchAddress("Jets_chargedMultiplicity", &Jets_chargedMultiplicity, &b_Jets_chargedMultiplicity);
   fChain->SetBranchAddress("Jets_electronEnergyFraction", &Jets_electronEnergyFraction, &b_Jets_electronEnergyFraction);
   fChain->SetBranchAddress("Jets_electronMultiplicity", &Jets_electronMultiplicity, &b_Jets_electronMultiplicity);
   fChain->SetBranchAddress("Jets_hadronFlavor", &Jets_hadronFlavor, &b_Jets_hadronFlavor);
   fChain->SetBranchAddress("Jets_hfEMEnergyFraction", &Jets_hfEMEnergyFraction, &b_Jets_hfEMEnergyFraction);
   fChain->SetBranchAddress("Jets_hfHadronEnergyFraction", &Jets_hfHadronEnergyFraction, &b_Jets_hfHadronEnergyFraction);
   fChain->SetBranchAddress("Jets_HTMask", &Jets_HTMask, &b_Jets_HTMask);
   fChain->SetBranchAddress("Jets_ID", &Jets_ID, &b_Jets_ID);
   fChain->SetBranchAddress("Jets_jecFactor", &Jets_jecFactor, &b_Jets_jecFactor);
   fChain->SetBranchAddress("Jets_jecUnc", &Jets_jecUnc, &b_Jets_jecUnc);
   fChain->SetBranchAddress("Jets_LeptonMask", &Jets_LeptonMask, &b_Jets_LeptonMask);
   fChain->SetBranchAddress("Jets_MHTMask", &Jets_MHTMask, &b_Jets_MHTMask);
   fChain->SetBranchAddress("Jets_multiplicity", &Jets_multiplicity, &b_Jets_multiplicity);
   fChain->SetBranchAddress("Jets_muonEnergyFraction", &Jets_muonEnergyFraction, &b_Jets_muonEnergyFraction);
   fChain->SetBranchAddress("Jets_muonMultiplicity", &Jets_muonMultiplicity, &b_Jets_muonMultiplicity);
   fChain->SetBranchAddress("Jets_neutralEmEnergyFraction", &Jets_neutralEmEnergyFraction, &b_Jets_neutralEmEnergyFraction);
   fChain->SetBranchAddress("Jets_neutralHadronEnergyFraction", &Jets_neutralHadronEnergyFraction, &b_Jets_neutralHadronEnergyFraction);
   fChain->SetBranchAddress("Jets_neutralHadronMultiplicity", &Jets_neutralHadronMultiplicity, &b_Jets_neutralHadronMultiplicity);
   fChain->SetBranchAddress("Jets_neutralMultiplicity", &Jets_neutralMultiplicity, &b_Jets_neutralMultiplicity);
   fChain->SetBranchAddress("Jets_partonFlavor", &Jets_partonFlavor, &b_Jets_partonFlavor);
   fChain->SetBranchAddress("Jets_photonEnergyFraction", &Jets_photonEnergyFraction, &b_Jets_photonEnergyFraction);
   fChain->SetBranchAddress("Jets_photonMultiplicity", &Jets_photonMultiplicity, &b_Jets_photonMultiplicity);
   fChain->SetBranchAddress("Jets_pileupId", &Jets_pileupId, &b_Jets_pileupId);
   fChain->SetBranchAddress("Jets_ptD", &Jets_ptD, &b_Jets_ptD);
   fChain->SetBranchAddress("Jets_qgLikelihood", &Jets_qgLikelihood, &b_Jets_qgLikelihood);
   fChain->SetBranchAddress("JetsAK15", &JetsAK15_, &b_JetsAK15_);
   fChain->SetBranchAddress("JetsAK15.fCoordinates.fPt", JetsAK15_fCoordinates_fPt, &b_JetsAK15_fCoordinates_fPt);
   fChain->SetBranchAddress("JetsAK15.fCoordinates.fEta", JetsAK15_fCoordinates_fEta, &b_JetsAK15_fCoordinates_fEta);
   fChain->SetBranchAddress("JetsAK15.fCoordinates.fPhi", JetsAK15_fCoordinates_fPhi, &b_JetsAK15_fCoordinates_fPhi);
   fChain->SetBranchAddress("JetsAK15.fCoordinates.fE", JetsAK15_fCoordinates_fE, &b_JetsAK15_fCoordinates_fE);
   fChain->SetBranchAddress("JetsAK15_axismajor", &JetsAK15_axismajor, &b_JetsAK15_axismajor);
   fChain->SetBranchAddress("JetsAK15_axisminor", &JetsAK15_axisminor, &b_JetsAK15_axisminor);
   fChain->SetBranchAddress("JetsAK15_chargedEmEnergyFraction", &JetsAK15_chargedEmEnergyFraction, &b_JetsAK15_chargedEmEnergyFraction);
   fChain->SetBranchAddress("JetsAK15_chargedHadronEnergyFraction", &JetsAK15_chargedHadronEnergyFraction, &b_JetsAK15_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("JetsAK15_chargedHadronMultiplicity", &JetsAK15_chargedHadronMultiplicity, &b_JetsAK15_chargedHadronMultiplicity);
   fChain->SetBranchAddress("JetsAK15_chargedMultiplicity", &JetsAK15_chargedMultiplicity, &b_JetsAK15_chargedMultiplicity);
   fChain->SetBranchAddress("JetsAK15_constituentsIndex", &JetsAK15_constituentsIndex, &b_JetsAK15_constituentsIndex);
   fChain->SetBranchAddress("JetsAK15_constituentsIndexCounts", &JetsAK15_constituentsIndexCounts, &b_JetsAK15_constituentsIndexCounts);
   fChain->SetBranchAddress("JetsAK15_DeepMassDecorrelTagbbvsLight", &JetsAK15_DeepMassDecorrelTagbbvsLight, &b_JetsAK15_DeepMassDecorrelTagbbvsLight);
   fChain->SetBranchAddress("JetsAK15_DeepMassDecorrelTagHbbvsQCD", &JetsAK15_DeepMassDecorrelTagHbbvsQCD, &b_JetsAK15_DeepMassDecorrelTagHbbvsQCD);
   fChain->SetBranchAddress("JetsAK15_DeepMassDecorrelTagTvsQCD", &JetsAK15_DeepMassDecorrelTagTvsQCD, &b_JetsAK15_DeepMassDecorrelTagTvsQCD);
   fChain->SetBranchAddress("JetsAK15_DeepMassDecorrelTagWvsQCD", &JetsAK15_DeepMassDecorrelTagWvsQCD, &b_JetsAK15_DeepMassDecorrelTagWvsQCD);
   fChain->SetBranchAddress("JetsAK15_DeepMassDecorrelTagZbbvsQCD", &JetsAK15_DeepMassDecorrelTagZbbvsQCD, &b_JetsAK15_DeepMassDecorrelTagZbbvsQCD);
   fChain->SetBranchAddress("JetsAK15_DeepMassDecorrelTagZHbbvsQCD", &JetsAK15_DeepMassDecorrelTagZHbbvsQCD, &b_JetsAK15_DeepMassDecorrelTagZHbbvsQCD);
   fChain->SetBranchAddress("JetsAK15_DeepMassDecorrelTagZvsQCD", &JetsAK15_DeepMassDecorrelTagZvsQCD, &b_JetsAK15_DeepMassDecorrelTagZvsQCD);
   fChain->SetBranchAddress("JetsAK15_DeepTagHbbvsQCD", &JetsAK15_DeepTagHbbvsQCD, &b_JetsAK15_DeepTagHbbvsQCD);
   fChain->SetBranchAddress("JetsAK15_DeepTagTvsQCD", &JetsAK15_DeepTagTvsQCD, &b_JetsAK15_DeepTagTvsQCD);
   fChain->SetBranchAddress("JetsAK15_DeepTagWvsQCD", &JetsAK15_DeepTagWvsQCD, &b_JetsAK15_DeepTagWvsQCD);
   fChain->SetBranchAddress("JetsAK15_DeepTagZbbvsQCD", &JetsAK15_DeepTagZbbvsQCD, &b_JetsAK15_DeepTagZbbvsQCD);
   fChain->SetBranchAddress("JetsAK15_DeepTagZvsQCD", &JetsAK15_DeepTagZvsQCD, &b_JetsAK15_DeepTagZvsQCD);
   fChain->SetBranchAddress("JetsAK15_doubleBDiscriminator", &JetsAK15_doubleBDiscriminator, &b_JetsAK15_doubleBDiscriminator);
   fChain->SetBranchAddress("JetsAK15_ecfC2b1", &JetsAK15_ecfC2b1, &b_JetsAK15_ecfC2b1);
   fChain->SetBranchAddress("JetsAK15_ecfC2b2", &JetsAK15_ecfC2b2, &b_JetsAK15_ecfC2b2);
   fChain->SetBranchAddress("JetsAK15_ecfD2b1", &JetsAK15_ecfD2b1, &b_JetsAK15_ecfD2b1);
   fChain->SetBranchAddress("JetsAK15_ecfD2b2", &JetsAK15_ecfD2b2, &b_JetsAK15_ecfD2b2);
   fChain->SetBranchAddress("JetsAK15_ecfM2b1", &JetsAK15_ecfM2b1, &b_JetsAK15_ecfM2b1);
   fChain->SetBranchAddress("JetsAK15_ecfM2b2", &JetsAK15_ecfM2b2, &b_JetsAK15_ecfM2b2);
   fChain->SetBranchAddress("JetsAK15_ecfN2b1", &JetsAK15_ecfN2b1, &b_JetsAK15_ecfN2b1);
   fChain->SetBranchAddress("JetsAK15_ecfN2b2", &JetsAK15_ecfN2b2, &b_JetsAK15_ecfN2b2);
   fChain->SetBranchAddress("JetsAK15_electronEnergyFraction", &JetsAK15_electronEnergyFraction, &b_JetsAK15_electronEnergyFraction);
   fChain->SetBranchAddress("JetsAK15_electronMultiplicity", &JetsAK15_electronMultiplicity, &b_JetsAK15_electronMultiplicity);
   fChain->SetBranchAddress("JetsAK15_girth", &JetsAK15_girth, &b_JetsAK15_girth);
   fChain->SetBranchAddress("JetsAK15_hadronFlavor", &JetsAK15_hadronFlavor, &b_JetsAK15_hadronFlavor);
   fChain->SetBranchAddress("JetsAK15_hfEMEnergyFraction", &JetsAK15_hfEMEnergyFraction, &b_JetsAK15_hfEMEnergyFraction);
   fChain->SetBranchAddress("JetsAK15_hfHadronEnergyFraction", &JetsAK15_hfHadronEnergyFraction, &b_JetsAK15_hfHadronEnergyFraction);
   fChain->SetBranchAddress("JetsAK15_ID", &JetsAK15_ID, &b_JetsAK15_ID);
   fChain->SetBranchAddress("JetsAK15_jecFactor", &JetsAK15_jecFactor, &b_JetsAK15_jecFactor);
   fChain->SetBranchAddress("JetsAK15_multiplicity", &JetsAK15_multiplicity, &b_JetsAK15_multiplicity);
   fChain->SetBranchAddress("JetsAK15_muonEnergyFraction", &JetsAK15_muonEnergyFraction, &b_JetsAK15_muonEnergyFraction);
   fChain->SetBranchAddress("JetsAK15_muonMultiplicity", &JetsAK15_muonMultiplicity, &b_JetsAK15_muonMultiplicity);
   fChain->SetBranchAddress("JetsAK15_neutralEmEnergyFraction", &JetsAK15_neutralEmEnergyFraction, &b_JetsAK15_neutralEmEnergyFraction);
   fChain->SetBranchAddress("JetsAK15_neutralHadronEnergyFraction", &JetsAK15_neutralHadronEnergyFraction, &b_JetsAK15_neutralHadronEnergyFraction);
   fChain->SetBranchAddress("JetsAK15_neutralHadronMultiplicity", &JetsAK15_neutralHadronMultiplicity, &b_JetsAK15_neutralHadronMultiplicity);
   fChain->SetBranchAddress("JetsAK15_neutralMultiplicity", &JetsAK15_neutralMultiplicity, &b_JetsAK15_neutralMultiplicity);
   fChain->SetBranchAddress("JetsAK15_NsubjettinessTau1", &JetsAK15_NsubjettinessTau1, &b_JetsAK15_NsubjettinessTau1);
   fChain->SetBranchAddress("JetsAK15_NsubjettinessTau2", &JetsAK15_NsubjettinessTau2, &b_JetsAK15_NsubjettinessTau2);
   fChain->SetBranchAddress("JetsAK15_NsubjettinessTau3", &JetsAK15_NsubjettinessTau3, &b_JetsAK15_NsubjettinessTau3);
   fChain->SetBranchAddress("JetsAK15_NsubjettinessTau4", &JetsAK15_NsubjettinessTau4, &b_JetsAK15_NsubjettinessTau4);
   fChain->SetBranchAddress("JetsAK15_NumBhadrons", &JetsAK15_NumBhadrons, &b_JetsAK15_NumBhadrons);
   fChain->SetBranchAddress("JetsAK15_NumChadrons", &JetsAK15_NumChadrons, &b_JetsAK15_NumChadrons);
   fChain->SetBranchAddress("JetsAK15_partonFlavor", &JetsAK15_partonFlavor, &b_JetsAK15_partonFlavor);
   fChain->SetBranchAddress("JetsAK15_pfMassIndependentDeepDoubleBvLJetTagsProbHbb", &JetsAK15_pfMassIndependentDeepDoubleBvLJetTagsProbHbb, &b_JetsAK15_pfMassIndependentDeepDoubleBvLJetTagsProbHbb);
   fChain->SetBranchAddress("JetsAK15_photonEnergyFraction", &JetsAK15_photonEnergyFraction, &b_JetsAK15_photonEnergyFraction);
   fChain->SetBranchAddress("JetsAK15_photonMultiplicity", &JetsAK15_photonMultiplicity, &b_JetsAK15_photonMultiplicity);
   fChain->SetBranchAddress("JetsAK15_ptD", &JetsAK15_ptD, &b_JetsAK15_ptD);
   fChain->SetBranchAddress("JetsAK15_softDropMass", &JetsAK15_softDropMass, &b_JetsAK15_softDropMass);
   fChain->SetBranchAddress("JetsAK15_softDropMassBeta1", &JetsAK15_softDropMassBeta1, &b_JetsAK15_softDropMassBeta1);
   fChain->SetBranchAddress("JetsAK15_subjets", &JetsAK15_subjets_, &b_JetsAK15_subjets_);
   fChain->SetBranchAddress("JetsAK15_subjets.fCoordinates.fPt", JetsAK15_subjets_fCoordinates_fPt, &b_JetsAK15_subjets_fCoordinates_fPt);
   fChain->SetBranchAddress("JetsAK15_subjets.fCoordinates.fEta", JetsAK15_subjets_fCoordinates_fEta, &b_JetsAK15_subjets_fCoordinates_fEta);
   fChain->SetBranchAddress("JetsAK15_subjets.fCoordinates.fPhi", JetsAK15_subjets_fCoordinates_fPhi, &b_JetsAK15_subjets_fCoordinates_fPhi);
   fChain->SetBranchAddress("JetsAK15_subjets.fCoordinates.fE", JetsAK15_subjets_fCoordinates_fE, &b_JetsAK15_subjets_fCoordinates_fE);
   fChain->SetBranchAddress("JetsAK15_subjetsCounts", &JetsAK15_subjetsCounts, &b_JetsAK15_subjetsCounts);
   fChain->SetBranchAddress("JetsAK8", &JetsAK8_, &b_JetsAK8_);
   fChain->SetBranchAddress("JetsAK8.fCoordinates.fPt", JetsAK8_fCoordinates_fPt, &b_JetsAK8_fCoordinates_fPt);
   fChain->SetBranchAddress("JetsAK8.fCoordinates.fEta", JetsAK8_fCoordinates_fEta, &b_JetsAK8_fCoordinates_fEta);
   fChain->SetBranchAddress("JetsAK8.fCoordinates.fPhi", JetsAK8_fCoordinates_fPhi, &b_JetsAK8_fCoordinates_fPhi);
   fChain->SetBranchAddress("JetsAK8.fCoordinates.fE", JetsAK8_fCoordinates_fE, &b_JetsAK8_fCoordinates_fE);
   fChain->SetBranchAddress("JetsAK8_axismajor", &JetsAK8_axismajor, &b_JetsAK8_axismajor);
   fChain->SetBranchAddress("JetsAK8_axisminor", &JetsAK8_axisminor, &b_JetsAK8_axisminor);
   fChain->SetBranchAddress("JetsAK8_chargedEmEnergyFraction", &JetsAK8_chargedEmEnergyFraction, &b_JetsAK8_chargedEmEnergyFraction);
   fChain->SetBranchAddress("JetsAK8_chargedHadronEnergyFraction", &JetsAK8_chargedHadronEnergyFraction, &b_JetsAK8_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("JetsAK8_chargedHadronMultiplicity", &JetsAK8_chargedHadronMultiplicity, &b_JetsAK8_chargedHadronMultiplicity);
   fChain->SetBranchAddress("JetsAK8_chargedMultiplicity", &JetsAK8_chargedMultiplicity, &b_JetsAK8_chargedMultiplicity);
   fChain->SetBranchAddress("JetsAK8_constituentsIndex", &JetsAK8_constituentsIndex, &b_JetsAK8_constituentsIndex);
   fChain->SetBranchAddress("JetsAK8_constituentsIndexCounts", &JetsAK8_constituentsIndexCounts, &b_JetsAK8_constituentsIndexCounts);
   fChain->SetBranchAddress("JetsAK8_DeepMassDecorrelTagbbvsLight", &JetsAK8_DeepMassDecorrelTagbbvsLight, &b_JetsAK8_DeepMassDecorrelTagbbvsLight);
   fChain->SetBranchAddress("JetsAK8_DeepMassDecorrelTagHbbvsQCD", &JetsAK8_DeepMassDecorrelTagHbbvsQCD, &b_JetsAK8_DeepMassDecorrelTagHbbvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepMassDecorrelTagTvsQCD", &JetsAK8_DeepMassDecorrelTagTvsQCD, &b_JetsAK8_DeepMassDecorrelTagTvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepMassDecorrelTagWvsQCD", &JetsAK8_DeepMassDecorrelTagWvsQCD, &b_JetsAK8_DeepMassDecorrelTagWvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepMassDecorrelTagZbbvsQCD", &JetsAK8_DeepMassDecorrelTagZbbvsQCD, &b_JetsAK8_DeepMassDecorrelTagZbbvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepMassDecorrelTagZHbbvsQCD", &JetsAK8_DeepMassDecorrelTagZHbbvsQCD, &b_JetsAK8_DeepMassDecorrelTagZHbbvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepMassDecorrelTagZvsQCD", &JetsAK8_DeepMassDecorrelTagZvsQCD, &b_JetsAK8_DeepMassDecorrelTagZvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepTagHbbvsQCD", &JetsAK8_DeepTagHbbvsQCD, &b_JetsAK8_DeepTagHbbvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepTagTvsQCD", &JetsAK8_DeepTagTvsQCD, &b_JetsAK8_DeepTagTvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepTagWvsQCD", &JetsAK8_DeepTagWvsQCD, &b_JetsAK8_DeepTagWvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepTagZbbvsQCD", &JetsAK8_DeepTagZbbvsQCD, &b_JetsAK8_DeepTagZbbvsQCD);
   fChain->SetBranchAddress("JetsAK8_DeepTagZvsQCD", &JetsAK8_DeepTagZvsQCD, &b_JetsAK8_DeepTagZvsQCD);
   fChain->SetBranchAddress("JetsAK8_doubleBDiscriminator", &JetsAK8_doubleBDiscriminator, &b_JetsAK8_doubleBDiscriminator);
   fChain->SetBranchAddress("JetsAK8_ecfN2b1", &JetsAK8_ecfN2b1, &b_JetsAK8_ecfN2b1);
   fChain->SetBranchAddress("JetsAK8_ecfN2b2", &JetsAK8_ecfN2b2, &b_JetsAK8_ecfN2b2);
   fChain->SetBranchAddress("JetsAK8_ecfN3b1", &JetsAK8_ecfN3b1, &b_JetsAK8_ecfN3b1);
   fChain->SetBranchAddress("JetsAK8_ecfN3b2", &JetsAK8_ecfN3b2, &b_JetsAK8_ecfN3b2);
   fChain->SetBranchAddress("JetsAK8_electronEnergyFraction", &JetsAK8_electronEnergyFraction, &b_JetsAK8_electronEnergyFraction);
   fChain->SetBranchAddress("JetsAK8_electronMultiplicity", &JetsAK8_electronMultiplicity, &b_JetsAK8_electronMultiplicity);
   fChain->SetBranchAddress("JetsAK8_girth", &JetsAK8_girth, &b_JetsAK8_girth);
   fChain->SetBranchAddress("JetsAK8_hadronFlavor", &JetsAK8_hadronFlavor, &b_JetsAK8_hadronFlavor);
   fChain->SetBranchAddress("JetsAK8_hfEMEnergyFraction", &JetsAK8_hfEMEnergyFraction, &b_JetsAK8_hfEMEnergyFraction);
   fChain->SetBranchAddress("JetsAK8_hfHadronEnergyFraction", &JetsAK8_hfHadronEnergyFraction, &b_JetsAK8_hfHadronEnergyFraction);
   fChain->SetBranchAddress("JetsAK8_ID", &JetsAK8_ID, &b_JetsAK8_ID);
   fChain->SetBranchAddress("JetsAK8_jecFactor", &JetsAK8_jecFactor, &b_JetsAK8_jecFactor);
   fChain->SetBranchAddress("JetsAK8_jecUnc", &JetsAK8_jecUnc, &b_JetsAK8_jecUnc);
   fChain->SetBranchAddress("JetsAK8_multiplicity", &JetsAK8_multiplicity, &b_JetsAK8_multiplicity);
   fChain->SetBranchAddress("JetsAK8_muonEnergyFraction", &JetsAK8_muonEnergyFraction, &b_JetsAK8_muonEnergyFraction);
   fChain->SetBranchAddress("JetsAK8_muonMultiplicity", &JetsAK8_muonMultiplicity, &b_JetsAK8_muonMultiplicity);
   fChain->SetBranchAddress("JetsAK8_neutralEmEnergyFraction", &JetsAK8_neutralEmEnergyFraction, &b_JetsAK8_neutralEmEnergyFraction);
   fChain->SetBranchAddress("JetsAK8_neutralHadronEnergyFraction", &JetsAK8_neutralHadronEnergyFraction, &b_JetsAK8_neutralHadronEnergyFraction);
   fChain->SetBranchAddress("JetsAK8_neutralHadronMultiplicity", &JetsAK8_neutralHadronMultiplicity, &b_JetsAK8_neutralHadronMultiplicity);
   fChain->SetBranchAddress("JetsAK8_neutralMultiplicity", &JetsAK8_neutralMultiplicity, &b_JetsAK8_neutralMultiplicity);
   fChain->SetBranchAddress("JetsAK8_NsubjettinessTau1", &JetsAK8_NsubjettinessTau1, &b_JetsAK8_NsubjettinessTau1);
   fChain->SetBranchAddress("JetsAK8_NsubjettinessTau2", &JetsAK8_NsubjettinessTau2, &b_JetsAK8_NsubjettinessTau2);
   fChain->SetBranchAddress("JetsAK8_NsubjettinessTau3", &JetsAK8_NsubjettinessTau3, &b_JetsAK8_NsubjettinessTau3);
   fChain->SetBranchAddress("JetsAK8_NumBhadrons", &JetsAK8_NumBhadrons, &b_JetsAK8_NumBhadrons);
   fChain->SetBranchAddress("JetsAK8_NumChadrons", &JetsAK8_NumChadrons, &b_JetsAK8_NumChadrons);
   fChain->SetBranchAddress("JetsAK8_partonFlavor", &JetsAK8_partonFlavor, &b_JetsAK8_partonFlavor);
   fChain->SetBranchAddress("JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb", &JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb, &b_JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb);
   fChain->SetBranchAddress("JetsAK8_photonEnergyFraction", &JetsAK8_photonEnergyFraction, &b_JetsAK8_photonEnergyFraction);
   fChain->SetBranchAddress("JetsAK8_photonMultiplicity", &JetsAK8_photonMultiplicity, &b_JetsAK8_photonMultiplicity);
   fChain->SetBranchAddress("JetsAK8_ptD", &JetsAK8_ptD, &b_JetsAK8_ptD);
   fChain->SetBranchAddress("JetsAK8_softDropMass", &JetsAK8_softDropMass, &b_JetsAK8_softDropMass);
   fChain->SetBranchAddress("JetsAK8_subjets", &JetsAK8_subjets_, &b_JetsAK8_subjets_);
   fChain->SetBranchAddress("JetsAK8_subjets.fCoordinates.fPt", JetsAK8_subjets_fCoordinates_fPt, &b_JetsAK8_subjets_fCoordinates_fPt);
   fChain->SetBranchAddress("JetsAK8_subjets.fCoordinates.fEta", JetsAK8_subjets_fCoordinates_fEta, &b_JetsAK8_subjets_fCoordinates_fEta);
   fChain->SetBranchAddress("JetsAK8_subjets.fCoordinates.fPhi", JetsAK8_subjets_fCoordinates_fPhi, &b_JetsAK8_subjets_fCoordinates_fPhi);
   fChain->SetBranchAddress("JetsAK8_subjets.fCoordinates.fE", JetsAK8_subjets_fCoordinates_fE, &b_JetsAK8_subjets_fCoordinates_fE);
   fChain->SetBranchAddress("JetsAK8_subjetsCounts", &JetsAK8_subjetsCounts, &b_JetsAK8_subjetsCounts);
   fChain->SetBranchAddress("JetsAK8_subjets_axismajor", &JetsAK8_subjets_axismajor, &b_JetsAK8_subjets_axismajor);
   fChain->SetBranchAddress("JetsAK8_subjets_axisminor", &JetsAK8_subjets_axisminor, &b_JetsAK8_subjets_axisminor);
   fChain->SetBranchAddress("JetsAK8_subjets_jecFactor", &JetsAK8_subjets_jecFactor, &b_JetsAK8_subjets_jecFactor);
   fChain->SetBranchAddress("JetsAK8_subjets_multiplicity", &JetsAK8_subjets_multiplicity, &b_JetsAK8_subjets_multiplicity);
   fChain->SetBranchAddress("JetsAK8_subjets_ptD", &JetsAK8_subjets_ptD, &b_JetsAK8_subjets_ptD);
   fChain->SetBranchAddress("JetsConstituents", &JetsConstituents_, &b_JetsConstituents_);
   fChain->SetBranchAddress("JetsConstituents.fCoordinates.fPt", JetsConstituents_fCoordinates_fPt, &b_JetsConstituents_fCoordinates_fPt);
   fChain->SetBranchAddress("JetsConstituents.fCoordinates.fEta", JetsConstituents_fCoordinates_fEta, &b_JetsConstituents_fCoordinates_fEta);
   fChain->SetBranchAddress("JetsConstituents.fCoordinates.fPhi", JetsConstituents_fCoordinates_fPhi, &b_JetsConstituents_fCoordinates_fPhi);
   fChain->SetBranchAddress("JetsConstituents.fCoordinates.fE", JetsConstituents_fCoordinates_fE, &b_JetsConstituents_fCoordinates_fE);
   fChain->SetBranchAddress("JetsConstituents_dxy", &JetsConstituents_dxy, &b_JetsConstituents_dxy);
   fChain->SetBranchAddress("JetsConstituents_dxysig", &JetsConstituents_dxysig, &b_JetsConstituents_dxysig);
   fChain->SetBranchAddress("JetsConstituents_dz", &JetsConstituents_dz, &b_JetsConstituents_dz);
   fChain->SetBranchAddress("JetsConstituents_dzsig", &JetsConstituents_dzsig, &b_JetsConstituents_dzsig);
   fChain->SetBranchAddress("JetsConstituents_PdgId", &JetsConstituents_PdgId, &b_JetsConstituents_PdgId);
   fChain->SetBranchAddress("JetsConstituents_PuppiWeight", &JetsConstituents_PuppiWeight, &b_JetsConstituents_PuppiWeight);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
   fChain->SetBranchAddress("METSignificance", &METSignificance, &b_METSignificance);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("MHTOnline", &MHTOnline, &b_MHTOnline);
   fChain->SetBranchAddress("MHTPhi", &MHTPhi, &b_MHTPhi);
   fChain->SetBranchAddress("MJJ_AK8", &MJJ_AK8, &b_MJJ_AK8);
   fChain->SetBranchAddress("Mmc_AK8", &Mmc_AK8, &b_Mmc_AK8);
   fChain->SetBranchAddress("MT_AK8", &MT_AK8, &b_MT_AK8);
   fChain->SetBranchAddress("Muons", &Muons_, &b_Muons_);
   fChain->SetBranchAddress("Muons.fCoordinates.fPt", Muons_fCoordinates_fPt, &b_Muons_fCoordinates_fPt);
   fChain->SetBranchAddress("Muons.fCoordinates.fEta", Muons_fCoordinates_fEta, &b_Muons_fCoordinates_fEta);
   fChain->SetBranchAddress("Muons.fCoordinates.fPhi", Muons_fCoordinates_fPhi, &b_Muons_fCoordinates_fPhi);
   fChain->SetBranchAddress("Muons.fCoordinates.fE", Muons_fCoordinates_fE, &b_Muons_fCoordinates_fE);
   fChain->SetBranchAddress("Muons_charge", &Muons_charge, &b_Muons_charge);
   fChain->SetBranchAddress("Muons_iso", &Muons_iso, &b_Muons_iso);
   fChain->SetBranchAddress("Muons_mediumID", &Muons_mediumID, &b_Muons_mediumID);
   fChain->SetBranchAddress("Muons_MTW", &Muons_MTW, &b_Muons_MTW);
   fChain->SetBranchAddress("Muons_passIso", &Muons_passIso, &b_Muons_passIso);
   fChain->SetBranchAddress("Muons_tightID", &Muons_tightID, &b_Muons_tightID);
   fChain->SetBranchAddress("nAllVertices", &nAllVertices, &b_nAllVertices);
   fChain->SetBranchAddress("NElectrons", &NElectrons, &b_NElectrons);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("NMuons", &NMuons, &b_NMuons);
   fChain->SetBranchAddress("NonPrefiringProb", &NonPrefiringProb, &b_NonPrefiringProb);
   fChain->SetBranchAddress("NonPrefiringProbDown", &NonPrefiringProbDown, &b_NonPrefiringProbDown);
   fChain->SetBranchAddress("NonPrefiringProbECAL", &NonPrefiringProbECAL, &b_NonPrefiringProbECAL);
   fChain->SetBranchAddress("NonPrefiringProbECALDown", &NonPrefiringProbECALDown, &b_NonPrefiringProbECALDown);
   fChain->SetBranchAddress("NonPrefiringProbECALUp", &NonPrefiringProbECALUp, &b_NonPrefiringProbECALUp);
   fChain->SetBranchAddress("NonPrefiringProbMuon", &NonPrefiringProbMuon, &b_NonPrefiringProbMuon);
   fChain->SetBranchAddress("NonPrefiringProbMuonDown", &NonPrefiringProbMuonDown, &b_NonPrefiringProbMuonDown);
   fChain->SetBranchAddress("NonPrefiringProbMuonUp", &NonPrefiringProbMuonUp, &b_NonPrefiringProbMuonUp);
   fChain->SetBranchAddress("NonPrefiringProbUp", &NonPrefiringProbUp, &b_NonPrefiringProbUp);
   fChain->SetBranchAddress("NVtx", &NVtx, &b_NVtx);
   fChain->SetBranchAddress("PFCaloMETRatio", &PFCaloMETRatio, &b_PFCaloMETRatio);
   fChain->SetBranchAddress("Photons", &Photons_, &b_Photons_);
   fChain->SetBranchAddress("Photons.fCoordinates.fPt", Photons_fCoordinates_fPt, &b_Photons_fCoordinates_fPt);
   fChain->SetBranchAddress("Photons.fCoordinates.fEta", Photons_fCoordinates_fEta, &b_Photons_fCoordinates_fEta);
   fChain->SetBranchAddress("Photons.fCoordinates.fPhi", Photons_fCoordinates_fPhi, &b_Photons_fCoordinates_fPhi);
   fChain->SetBranchAddress("Photons.fCoordinates.fE", Photons_fCoordinates_fE, &b_Photons_fCoordinates_fE);
   fChain->SetBranchAddress("Photons_electronFakes", &Photons_electronFakes, &b_Photons_electronFakes);
   fChain->SetBranchAddress("Photons_fullID", &Photons_fullID, &b_Photons_fullID);
   fChain->SetBranchAddress("Photons_genMatched", &Photons_genMatched, &b_Photons_genMatched);
   fChain->SetBranchAddress("Photons_hadTowOverEM", &Photons_hadTowOverEM, &b_Photons_hadTowOverEM);
   fChain->SetBranchAddress("Photons_hasPixelSeed", &Photons_hasPixelSeed, &b_Photons_hasPixelSeed);
   fChain->SetBranchAddress("Photons_isEB", &Photons_isEB, &b_Photons_isEB);
   fChain->SetBranchAddress("Photons_nonPrompt", &Photons_nonPrompt, &b_Photons_nonPrompt);
   fChain->SetBranchAddress("Photons_passElectronVeto", &Photons_passElectronVeto, &b_Photons_passElectronVeto);
   fChain->SetBranchAddress("Photons_pfChargedIso", &Photons_pfChargedIso, &b_Photons_pfChargedIso);
   fChain->SetBranchAddress("Photons_pfChargedIsoRhoCorr", &Photons_pfChargedIsoRhoCorr, &b_Photons_pfChargedIsoRhoCorr);
   fChain->SetBranchAddress("Photons_pfGammaIso", &Photons_pfGammaIso, &b_Photons_pfGammaIso);
   fChain->SetBranchAddress("Photons_pfGammaIsoRhoCorr", &Photons_pfGammaIsoRhoCorr, &b_Photons_pfGammaIsoRhoCorr);
   fChain->SetBranchAddress("Photons_pfNeutralIso", &Photons_pfNeutralIso, &b_Photons_pfNeutralIso);
   fChain->SetBranchAddress("Photons_pfNeutralIsoRhoCorr", &Photons_pfNeutralIsoRhoCorr, &b_Photons_pfNeutralIsoRhoCorr);
   fChain->SetBranchAddress("Photons_sigmaIetaIeta", &Photons_sigmaIetaIeta, &b_Photons_sigmaIetaIeta);
   fChain->SetBranchAddress("PrescaleWeightHT", &PrescaleWeightHT, &b_PrescaleWeightHT);
   fChain->SetBranchAddress("PrimaryVertexFilter", &PrimaryVertexFilter, &b_PrimaryVertexFilter);
   fChain->SetBranchAddress("TAPElectronTracks", &TAPElectronTracks_, &b_TAPElectronTracks_);
   fChain->SetBranchAddress("TAPElectronTracks.fCoordinates.fPt", TAPElectronTracks_fCoordinates_fPt, &b_TAPElectronTracks_fCoordinates_fPt);
   fChain->SetBranchAddress("TAPElectronTracks.fCoordinates.fEta", TAPElectronTracks_fCoordinates_fEta, &b_TAPElectronTracks_fCoordinates_fEta);
   fChain->SetBranchAddress("TAPElectronTracks.fCoordinates.fPhi", TAPElectronTracks_fCoordinates_fPhi, &b_TAPElectronTracks_fCoordinates_fPhi);
   fChain->SetBranchAddress("TAPElectronTracks.fCoordinates.fE", TAPElectronTracks_fCoordinates_fE, &b_TAPElectronTracks_fCoordinates_fE);
   fChain->SetBranchAddress("TAPElectronTracks_dxypv", &TAPElectronTracks_dxypv, &b_TAPElectronTracks_dxypv);
   fChain->SetBranchAddress("TAPElectronTracks_leptonMatch", &TAPElectronTracks_leptonMatch, &b_TAPElectronTracks_leptonMatch);
   fChain->SetBranchAddress("TAPElectronTracks_mT", &TAPElectronTracks_mT, &b_TAPElectronTracks_mT);
   fChain->SetBranchAddress("TAPElectronTracks_pfRelIso03chg", &TAPElectronTracks_pfRelIso03chg, &b_TAPElectronTracks_pfRelIso03chg);
   fChain->SetBranchAddress("TAPElectronTracks_trkiso", &TAPElectronTracks_trkiso, &b_TAPElectronTracks_trkiso);
   fChain->SetBranchAddress("TAPMuonTracks", &TAPMuonTracks_, &b_TAPMuonTracks_);
   fChain->SetBranchAddress("TAPMuonTracks.fCoordinates.fPt", TAPMuonTracks_fCoordinates_fPt, &b_TAPMuonTracks_fCoordinates_fPt);
   fChain->SetBranchAddress("TAPMuonTracks.fCoordinates.fEta", TAPMuonTracks_fCoordinates_fEta, &b_TAPMuonTracks_fCoordinates_fEta);
   fChain->SetBranchAddress("TAPMuonTracks.fCoordinates.fPhi", TAPMuonTracks_fCoordinates_fPhi, &b_TAPMuonTracks_fCoordinates_fPhi);
   fChain->SetBranchAddress("TAPMuonTracks.fCoordinates.fE", TAPMuonTracks_fCoordinates_fE, &b_TAPMuonTracks_fCoordinates_fE);
   fChain->SetBranchAddress("TAPMuonTracks_dxypv", &TAPMuonTracks_dxypv, &b_TAPMuonTracks_dxypv);
   fChain->SetBranchAddress("TAPMuonTracks_leptonMatch", &TAPMuonTracks_leptonMatch, &b_TAPMuonTracks_leptonMatch);
   fChain->SetBranchAddress("TAPMuonTracks_mT", &TAPMuonTracks_mT, &b_TAPMuonTracks_mT);
   fChain->SetBranchAddress("TAPMuonTracks_pfRelIso03chg", &TAPMuonTracks_pfRelIso03chg, &b_TAPMuonTracks_pfRelIso03chg);
   fChain->SetBranchAddress("TAPMuonTracks_trkiso", &TAPMuonTracks_trkiso, &b_TAPMuonTracks_trkiso);
   fChain->SetBranchAddress("TAPPionTracks", &TAPPionTracks_, &b_TAPPionTracks_);
   fChain->SetBranchAddress("TAPPionTracks.fCoordinates.fPt", TAPPionTracks_fCoordinates_fPt, &b_TAPPionTracks_fCoordinates_fPt);
   fChain->SetBranchAddress("TAPPionTracks.fCoordinates.fEta", TAPPionTracks_fCoordinates_fEta, &b_TAPPionTracks_fCoordinates_fEta);
   fChain->SetBranchAddress("TAPPionTracks.fCoordinates.fPhi", TAPPionTracks_fCoordinates_fPhi, &b_TAPPionTracks_fCoordinates_fPhi);
   fChain->SetBranchAddress("TAPPionTracks.fCoordinates.fE", TAPPionTracks_fCoordinates_fE, &b_TAPPionTracks_fCoordinates_fE);
   fChain->SetBranchAddress("TAPPionTracks_dxypv", &TAPPionTracks_dxypv, &b_TAPPionTracks_dxypv);
   fChain->SetBranchAddress("TAPPionTracks_leptonMatch", &TAPPionTracks_leptonMatch, &b_TAPPionTracks_leptonMatch);
   fChain->SetBranchAddress("TAPPionTracks_mT", &TAPPionTracks_mT, &b_TAPPionTracks_mT);
   fChain->SetBranchAddress("TAPPionTracks_pfRelIso03chg", &TAPPionTracks_pfRelIso03chg, &b_TAPPionTracks_pfRelIso03chg);
   fChain->SetBranchAddress("TAPPionTracks_trkiso", &TAPPionTracks_trkiso, &b_TAPPionTracks_trkiso);
   fChain->SetBranchAddress("TriggerPass", &TriggerPass, &b_TriggerPass);
   fChain->SetBranchAddress("TriggerPrescales", &TriggerPrescales, &b_TriggerPrescales);
   fChain->SetBranchAddress("TriggerVersion", &TriggerVersion, &b_TriggerVersion);
   Notify();
}

Bool_t allin1::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void allin1::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t allin1::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef allin1_cxx
