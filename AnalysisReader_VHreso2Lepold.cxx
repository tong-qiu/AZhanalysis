#include <type_traits>

#include "CxAODReader/EasyTree.h"
#include "CxAODReader_VHbb/AnalysisReader_VHQQ2Lep.h"
#include "CxAODReader_VHbb/MVATree_BoostedVHbb.h"
#include "CxAODReader_VHbb/MVATree_VHbb.h"
#include "CxAODTools_VHbb/TriggerTool_VHbb2lep.h"
#include "CxAODTools_VHbb/VHbb2lepEvtSelection.h"
#include "KinematicFit/KinematicFit.h"

#include <CxAODReader_VHbb/AnalysisReader_VHreso2Lep.h>
#include <TMVA/Reader.h>
#include "CxAODReader_VHbb/AnalysisReader_VHQQ.h"
#include "TXMLEngine.h"

#include <iomanip>

#define length(array) (sizeof(array) / sizeof(*(array)))

// default constructor for all the event variables
// allows to easily reset all values at the beginning of every run_2Lep_analysis
// call
AnalysisReader_VHreso2Lep::EventVariables::EventVariables()
    : el1(nullptr),
      el2(nullptr),
      mu1(nullptr),
      mu2(nullptr),
      met(nullptr),
      signalJets(),
      forwardJets(),
      fatJets(),
      trackJets(),
      taus(),
      trackJetsInLeadFJ(),
      trackJetsNotInLeadFJ(),
      selectedJets(),
      trackJetsForBTagging(),
      bTaggedTrackJetsNotInLeadFJ(),
      trackjetsForLabelling(),
      j1(),
      j2(),
      j3(),
      fj1(),
      fj2(),
      fj3(),
      nTaggedTrkJets(0),
      ZVec(),
      bbjVec(),
      HVec(),
      HVec_gsc(),
      HVec_onemu(),
      HVecPtReco(),
      HVec_ptreco(),
      HVec_kf(),
      VHVec(),
      VHVec_gsc(),
      VHVec_onemu(),
      VHVec_ptreco(),
      VHVec_kf(),
      lep1(),
      lep2(),
      resolvedH(),
      mergedH(),
      resolvedVH(),
      mergedVH(),
      HTsmallR(0),
      HTlargeR(0),
      METHTsmallR(0),
      METHTlargeR(0),
      eventFlag(0),
      btagWeight(1),
      triggerSF(1),
      JVTweight(1) {}

AnalysisReader_VHreso2Lep::AnalysisReader_VHreso2Lep()
    : AnalysisReader_VHQQ2Lep(), m_writeEasyTree(false), m_vars() {
  m_doMergeModel = false;
  m_truthLabeling = "TrackJetHybrid";
  m_bdtFileDir = "$WorkDir_DIR/data/CxAODReader_VHbb/BDT_2lep/v2";
  m_doKF_allcat = false;
}

AnalysisReader_VHreso2Lep::~AnalysisReader_VHreso2Lep() {}

EL::StatusCode AnalysisReader_VHreso2Lep::run_2Lep_analysis() {
  m_config->getif<bool>("writeEasyTree", m_writeEasyTree);

  // remove overlap between dilep and inclusive / non-all had ttbar and single
  // top Wt samples
  //----------------------------------------------------------------------------
  // requires to include both in the sample list!
  // apply to which samples? just to the nominal PP8 ones?
  // todo: move vetoDilepTtbarEvents implementation to VHQQ so that it can also
  // be used by VHcc? - wait until harmonised treatment available for ttbar and
  // single top?
  if (m_doRemoveDilepOverlap) {
    // ttbar
    if (m_mcChannel == 410470) {  // 410470: PP8 non-all-had ttbar
      int veto = vetoDilepTtbarEvents();
      if (veto == -1) {
        Error("run_2Lep_analysis()", "Props::codeTTBarDecay doesn't exist!");
        return EL::StatusCode::FAILURE;
      } else if (veto == 1) {
        // Event has to be vetoed
        return EL::StatusCode::SUCCESS;
      }
    }

    // single top Wt
    if (m_mcChannel == 410646 ||
        m_mcChannel ==
            410647) {  // 410646(7): PP8 inclusive single (anti-)top Wt
      int veto = vetoDilepWtEvents();
      // Event has to be vetoed
      if (veto == 1) return EL::StatusCode::SUCCESS;
    }
  }
  //----------------------------

  // reset physics metadata
  m_physicsMeta = PhysicsMetadata();

  // turn truth tagging off for merged analysis for the moment
  if (m_analysisStrategy != "Resolved") m_doTruthTagging = false;

  // EasyTree reset
  m_etree->Reset();
  m_etree->SetVariation(m_currentVar);
  m_vars = {};

  ResultVHbb2lep selectionResult =
      ((VHbb2lepEvtSelection *)m_eventSelection)->result();

  m_physicsMeta.channel = PhysicsMetadata::Channel::TwoLep;

  m_vars.el1 = selectionResult.el1;
  m_vars.el2 = selectionResult.el2;
  m_vars.mu1 = selectionResult.mu1;
  m_vars.mu2 = selectionResult.mu2;
  m_vars.met = selectionResult.met;
  m_vars.signalJets = selectionResult.signalJets;
  m_vars.forwardJets = selectionResult.forwardJets;
  m_vars.fatJets = selectionResult.fatJets;
  m_vars.trackJets = selectionResult.trackJets;
  m_vars.taus = selectionResult.taus;
  setevent_nJets(m_vars.signalJets, m_vars.forwardJets);

  ////////////////////////////////////////////////////////////
  //         OBJECT DEFINITIONS                             //
  ////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // **Leptons Definition** :
  EL_CHECK("run_2Lep_analysis()", setLeptonVariables());
  EL_CHECK("run_2Lep_analysis()", setTwoLeptonFlavour());

  bool isMu = (m_physicsMeta.flavor == PhysicsMetadata::Flavor::MuMu);
  bool isE = (m_physicsMeta.flavor == PhysicsMetadata::Flavor::ElEl);

  /////////////////////////////////////////////////////////////
  // **Definition** : Z
  m_vars.ZVec = m_vars.lep1.vec + m_vars.lep2.vec;

  /////////////////////////////////////////////////////////////
  // **Definition** : b-tagging

  // TODO: Harmonize large-R jet / track-jet association with resonant analysis
  if (m_doTruthTagging) {
    compute_TRF_tagging(m_vars.signalJets);
  } else {
    compute_btagging();
  }
  compute_fatjetTags(m_vars.trackJets, m_vars.fatJets,
                     &m_vars.trackJetsInLeadFJ, &m_vars.trackJetsNotInLeadFJ);

  /////////////////////////////////////////////////////////////
  // **Definition** : selected jets
  int tagcatExcl = -1;

  tagjet_selection(m_vars.signalJets, m_vars.forwardJets, m_vars.selectedJets,
                   tagcatExcl);
  if (m_model == Model::CUT || m_model == Model::MVA) {
    // if Forward jet has higher pT assign it as a third jet.
    if (m_vars.selectedJets.size() == 3) {
      for (unsigned int ifjet = 0; ifjet < m_vars.forwardJets.size(); ++ifjet) {
        if (m_vars.selectedJets.at(2)->pt() <
            m_vars.forwardJets.at(ifjet)->pt())
          m_vars.selectedJets.at(2) = m_vars.forwardJets.at(ifjet);
      }
    }
  }
  m_physicsMeta.nTags = tagcatExcl;
  // btags inside/outside fat jet
  if (m_vars.fatJets.size() > 0) {
    // number of lead / sublead b-tagged track jets associated to leading fat
    // jet
    m_physicsMeta.nTagsInFJ = Props::nBTags.get(m_vars.fatJets.at(0));
    ;
    // number of b-tagged track jets not associated to leading fat jet
    m_physicsMeta.nAddBTrkJets = Props::nAddBTags.get(m_vars.fatJets.at(0));
  }
  // number of b-tagged track jets
  m_vars.nTaggedTrkJets = 0;
  for (auto jet : m_vars.trackJets) {
    if (BTagProps::isTagged.get(jet)) m_vars.nTaggedTrkJets++;
  }

  // define jets for b-tagging weight computation and event labelling
  //------------------------------------------------------------------
  // consider two leading track jets matched to the lead. fat jet
  // and all track jets not matched to fat jet in b-tagging event weight
  // computation
  m_vars.trackJetsForBTagging = m_vars.trackJetsInLeadFJ;
  m_vars.trackJetsForBTagging.insert(m_vars.trackJetsForBTagging.end(),
                                     m_vars.trackJetsNotInLeadFJ.begin(),
                                     m_vars.trackJetsNotInLeadFJ.end());
  // find track jets used for flavour labelling the event: one (two) leading
  // ones inside the fat jet, leading (b-tagged) one outside of the fat jet
  // first: find b-tagged unmatched track jets
  for (auto jet : m_vars.trackJetsNotInLeadFJ) {
    if (BTagProps::isTagged.get(jet))
      m_vars.bTaggedTrackJetsNotInLeadFJ.push_back(jet);
  }
  // second: check consistency with what was derived before
  if ((unsigned int)m_physicsMeta.nAddBTrkJets !=
      m_vars.bTaggedTrackJetsNotInLeadFJ.size()) {
    Error("AnalysisReader_VHreso2Lep::run_2Lep_analysis",
          "m_physicsMeta.nAddBTrkJets != "
          "m_vars.bTaggedTrackJetsNotInLeadFJ.size()");
    return EL::StatusCode::FAILURE;
  }
  // third: select jets for event flavour labelling
  m_vars.trackjetsForLabelling = m_vars.trackJetsInLeadFJ;
  if (m_physicsMeta.nAddBTrkJets)
    m_vars.trackjetsForLabelling.insert(
        m_vars.trackjetsForLabelling.end(),
        m_vars.bTaggedTrackJetsNotInLeadFJ.begin(),
        m_vars.bTaggedTrackJetsNotInLeadFJ.end());
  else
    m_vars.trackjetsForLabelling.insert(m_vars.trackjetsForLabelling.end(),
                                        m_vars.trackJetsNotInLeadFJ.begin(),
                                        m_vars.trackJetsNotInLeadFJ.end());

  // define jets to be used for mBB + possible additional jet
  //---------------------------------------------------------
  EL_CHECK(
      "AnalysisReader_VHreso2Lep::run_2Lep_analysis()",
      setJetVariables(m_vars.j1, m_vars.j2, m_vars.j3, m_vars.selectedJets));
  EL_CHECK(
      "AnalysisReader_VHreso2Lep::run_2Lep_analysis()",
      setFatJetVariables(m_vars.fj1, m_vars.j2, m_vars.j3, m_vars.fatJets));

  /////////////////////////////////////////////////////////////
  // **Definition**: MET
  /*
   * unused as of 2018-12-20, but might be included in full Run2 paper
  double metSig = -1.;
  double metSig_PU = -1.;
  double metSig_soft = -1.;
  double metSig_hard = -1.;

  metSig = Props::metSig.get(m_vars.met);
  if (Props::metSig_PU.exists(m_vars.met)) {
    metSig_PU = Props::metSig_PU.get(m_vars.met);
  } else {
    metSig_hard = Props::metSig_hard.get(m_vars.met);
    metSig_soft = Props::metSig_soft.get(m_vars.met);
  }
  */
  /////////////////////////////////////////////////////////////
  // **Definition** : HT and METHT
  computeHT(m_vars.HTsmallR, m_vars.HTlargeR, m_vars.lep1, m_vars.lep2,
            m_vars.signalJets, m_vars.forwardJets, m_vars.fatJets);

  m_vars.METHTsmallR = m_vars.met->met() / sqrt(m_vars.HTsmallR);
  m_vars.METHTlargeR = m_vars.met->met() / sqrt(m_vars.HTlargeR);

  /////////////////////////////////////////////////////////////
  // **Definition** : m_vars.HVec
  setHiggsCandidate(m_vars.resolvedH, m_vars.mergedH, m_vars.j1, m_vars.j2,
                    m_vars.j1);
  bool kf_converged;
  if (m_doKF_allcat) {  // for debugging
    EL_CHECK("AnalysisReader_VHreso2Lep::run_2Lep_analysis()",
             getKFResolved(m_vars.selectedJets.size(), m_vars.j1, m_vars.j2,
                           m_vars.j3, m_vars.lep1, m_vars.lep2,
                           kf_converged));  // Call KF even when it is not used
                                            // for SM VHbb analysis
  }
  if ((m_physicsMeta.flavor == PhysicsMetadata::Flavor::ElEl ||
       m_physicsMeta.flavor == PhysicsMetadata::Flavor::MuMu ||
       m_physicsMeta.flavor == PhysicsMetadata::Flavor::ElMu ||
       m_physicsMeta.flavor == PhysicsMetadata::Flavor::MuEl) &&
      ((m_physicsMeta.nJets == 2) || (m_physicsMeta.nJets == 3)) &&
      m_physicsMeta.nSigJet >= 2 && m_physicsMeta.nTags == 2
      //&& kf_converged //Will be switched on for rel21 analysis
  ) {                      // if regions for KF
    if (!m_doKF_allcat) {  // KF has not been called yet
      EL_CHECK(
          "AnalysisReader_VHreso2Lep::run_2Lep_analysis()",
          getKFResolved(m_vars.selectedJets.size(), m_vars.j1, m_vars.j2,
                        m_vars.j3, m_vars.lep1, m_vars.lep2, kf_converged));
    }
    if (m_model == Model::CUT || m_model == Model::MVA) {
      m_vars.j1.vec_corr = m_vars.j1.vec_kf;
      m_vars.j2.vec_corr = m_vars.j2.vec_kf;
    }
  }
  m_vars.resolvedH.vec_kf = m_vars.j1.vec_kf + m_vars.j2.vec_kf;

  if ((m_physicsMeta.flavor == PhysicsMetadata::Flavor::ElEl ||
       m_physicsMeta.flavor == PhysicsMetadata::Flavor::MuMu ||
       m_physicsMeta.flavor == PhysicsMetadata::Flavor::ElMu ||
       m_physicsMeta.flavor == PhysicsMetadata::Flavor::MuEl) &&
      (m_vars.fatJets.size() == 1 &&
       ((m_physicsMeta.nJets == 0 && m_vars.selectedJets.size() == 0) ||
        (m_physicsMeta.nJets == 1 && m_vars.selectedJets.size() == 1) ||
        (m_physicsMeta.nJets == 2 && m_vars.selectedJets.size() == 2) ||
        (m_physicsMeta.nJets == 3 && m_vars.selectedJets.size() == 3)))) {
    EL_CHECK("AnalysisReader_VHreso2Lep::run_2Lep_analysis()",
             getKFMerged(m_physicsMeta.nJets, m_vars.j1, m_vars.lep1,
                         m_vars.lep2, m_vars.j1, m_vars.j2, m_vars.j3));
  }
  m_vars.mergedH.vec_kf = m_vars.j1.vec_kf;

  if (m_model == Model::CUT || m_model == Model::MVA) {
    setHiggsCandidate(m_vars.resolvedH, m_vars.mergedH, m_vars.j1, m_vars.j2,
                      m_vars.j1);
  }

  /////////////////////////////////////////////////////////////
  //  **Definition** : ZH system
  setVHCandidate(m_vars.resolvedVH, m_vars.mergedVH, m_vars.lep1, m_vars.lep2,
                 m_vars.j1, m_vars.j2, m_vars.fj1);

  ////////////////////////////////////////////////////////////
  //         CHECK SELECTION - STORE INFORMATION            //
  ////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////
  // C0: All events (CxAOD)
  //------------------------
  updateFlag(m_vars.eventFlag, DiLeptonCuts::AllCxAOD);

  // **Selection & SF**: trigger
  //----------------------------
  // C1: trigger
  if (pass2LepTrigger(m_vars.triggerSF, selectionResult))
    updateFlag(m_vars.eventFlag, DiLeptonCuts::Trigger);
  if (m_isMC) m_weight *= m_vars.triggerSF;

  // **Selection** : C2: same flavor
  //----------------------------
  if (isMu || isE)
    updateFlag(m_vars.eventFlag, DiLeptonCuts::SFLeptons);
  else
    updateFlag(m_vars.eventFlag, DiLeptonCuts::DFLeptons);

  if (m_isMC) m_weight *= Props::leptonSF.get(m_eventInfo);

  // **CR Selection** : cut on lepton pT
  //------------------------------------
  if (m_vars.lep1.vec.Pt() / 1000. >= 27. && m_vars.lep2.vec.Pt() / 1000. >= 20.) {
    updateFlag(m_vars.eventFlag, DiLeptonCuts::LeadLepPt27);
    if (m_vars.lep2.vec.Pt() / 1000. >= 25.)
      updateFlag(m_vars.eventFlag, DiLeptonCuts::pTBothLeptons25);
  }

  // **CR Selection** : cut on lepton charge and muon eta (for merged analysis)
  //-----------------------------------------------------
  if (isE) {
    updateFlag(m_vars.eventFlag, DiLeptonCuts::OSLeptons);
    updateFlag(m_vars.eventFlag, DiLeptonCuts::MuonEtaLT2p5);
  } else if (isMu) {
    if (m_vars.lep1.charge != m_vars.lep2.charge)
      updateFlag(m_vars.eventFlag, DiLeptonCuts::OSLeptons);
    if (fabs(m_vars.lep1.vec.Eta()) < 2.5 && fabs(m_vars.lep2.vec.Eta()) < 2.5)
      updateFlag(m_vars.eventFlag, DiLeptonCuts::MuonEtaLT2p5);
  } else {  // emu channel
    if (m_vars.lep1.charge != m_vars.lep2.charge)
      updateFlag(m_vars.eventFlag, DiLeptonCuts::OSLeptons);
    if (((m_vars.lep1.flav == lepFlav::mu) &&
         (fabs(m_vars.lep1.vec.Eta()) < 2.5)) ||
        ((m_vars.lep2.flav == lepFlav::mu) &&
         (fabs(m_vars.lep2.vec.Eta()) < 2.5)))
      updateFlag(m_vars.eventFlag, DiLeptonCuts::MuonEtaLT2p5);
  }

  // **Selection** : tau veto
  //----------------------------------------
  if (!m_vars.taus.size()) updateFlag(m_vars.eventFlag, DiLeptonCuts::tauVeto);

  // **Selection** : C3: Mll cut
  //-----------------------------
  bool passMllResolved(false), passMllMerged(false);
  if (m_model == Model::AZh || m_model == Model::HVT)
    checkMllCut(m_vars.ZVec, m_vars.resolvedVH.vec_corr,
                m_vars.mergedVH.vec_corr, passMllResolved, passMllMerged);
  else if (m_model == Model::CUT || m_model == Model::MVA) {
    if ((m_vars.ZVec.M() / 1000. >= 81.) && (m_vars.ZVec.M() / 1000. < 101.)) {
      passMllResolved = true;
      passMllMerged = true;
    }
  }

  if (passMllResolved)
    updateFlag(m_vars.eventFlag, DiLeptonCuts::MllZwindowResolved);
  if (passMllMerged)
    updateFlag(m_vars.eventFlag, DiLeptonCuts::MllZwindowMerged);
  if (m_vars.ZVec.M() / 1.e3 > 40.)
    updateFlag(m_vars.eventFlag, DiLeptonCuts::MllAbove40);

  //   **Selection** : C4: METHT cut
  //----------------------------------
  if (m_model == Model::AZh || m_model == Model::HVT) {
    if (m_vars.METHTsmallR / sqrt(1000.) <
        METHTcut(m_vars.resolvedVH.vec_corr.M()))
      updateFlag(m_vars.eventFlag, DiLeptonCuts::METHTResolved);
    if (m_vars.METHTsmallR / sqrt(1000.) <
        METHTcut(m_vars.mergedVH.vec_corr.M()))
      updateFlag(m_vars.eventFlag, DiLeptonCuts::METHTMerged);
  } else if (m_model == Model::CUT || m_model == Model::MVA) {
    if (m_vars.METHTsmallR / sqrt(1000.) < 3.5)
      updateFlag(m_vars.eventFlag, DiLeptonCuts::METHTResolved);
    updateFlag(m_vars.eventFlag, DiLeptonCuts::METHTMerged);
  }

  // **Selection** : C5: cut on nJet
  //-----------------------------------
  if (m_physicsMeta.nJets >= 2)
    updateFlag(m_vars.eventFlag, DiLeptonCuts::AtLeast2Jets);

  // **Selection** : C6: cut on nSignalJet
  //-----------------------------------------
  if (m_physicsMeta.nSigJet >= 2)
    updateFlag(m_vars.eventFlag, DiLeptonCuts::AtLeast2SigJets);
  if (m_isMC) {
    m_vars.JVTweight = compute_JVTSF(m_vars.signalJets);
    m_weight *= m_vars.JVTweight;
  }

  // **Selection** : B-tagging/Jet selection
  //-----------------------------------------
  // C7: at least 1 b
  if (m_physicsMeta.nTags >= 1)
    updateFlag(m_vars.eventFlag, DiLeptonCuts::AtLeast1B);
  // C8: exactly 2 b
  if (m_physicsMeta.nTags == 2)
    updateFlag(m_vars.eventFlag, DiLeptonCuts::Exactly2B);
  // if (m_physicsMeta.nTags < 3) updateFlag(m_vars.eventFlag,
  // DiLeptonCuts::Veto3B);

  // **Selection** :  C9: pt(leading jet) >= 45 GeV
  //-----------------------------------------
  if (m_vars.selectedJets.size() > 0 && m_vars.j1.vec.Pt() / 1.e3 > 45.)
    updateFlag(m_vars.eventFlag, DiLeptonCuts::PtB145);

  // **Selection** : dRCut for SM CBA
  //-----------------------------------------
  float dRCut = 3.0;
  if (m_vars.ZVec.Pt() * 0.001 > 200) {
    dRCut = 1.2;
  } else if (m_vars.ZVec.Pt() * 0.001 > 150) {
    dRCut = 1.8;
  }
  if (m_model == Model::CUT && m_vars.selectedJets.size() > 1 &&
      m_vars.j1.vec_onemu.DeltaR(m_vars.j2.vec_onemu) < dRCut) {
    updateFlag(m_vars.eventFlag, DiLeptonCuts::dRCut);
  }

  //  **Selection** : C10': mbb w/o mu-in-jet (+ pT reco) corr
  //------------------------------------------
  bool passMHResolved(false), passMHMerged(false), inMHSideBandResolved(false),
      inMHSideBandMerged(false);
  checkMbbWindow(m_vars.resolvedH.vec, m_vars.mergedH.vec, passMHResolved,
                 passMHMerged, inMHSideBandResolved, inMHSideBandMerged);

  if (passMHResolved)
    updateFlag(m_vars.eventFlag, DiLeptonCuts::mHNoCorrResolved);
  if (passMHMerged) updateFlag(m_vars.eventFlag, DiLeptonCuts::mHNoCorrMerged);

  // C10: mbb w/ mu-in-jet (+ pT reco) corr
  //----------------------------------------
  checkMbbWindow(m_vars.resolvedH.vec_corr, m_vars.mergedH.vec_corr,
                 passMHResolved, passMHMerged, inMHSideBandResolved,
                 inMHSideBandMerged);

  if (passMHResolved)
    updateFlag(m_vars.eventFlag, DiLeptonCuts::mHCorrResolved);
  if (passMHMerged) updateFlag(m_vars.eventFlag, DiLeptonCuts::mHCorrMerged);
  if (inMHSideBandResolved)
    updateFlag(m_vars.eventFlag, DiLeptonCuts::mHCorrInsideSideBandResolved);
  if (inMHSideBandMerged)
    updateFlag(m_vars.eventFlag, DiLeptonCuts::mHCorrInsideSideBandMerged);

  // **Selection** :  PtZ
  //------------------------
  if (m_vars.ZVec.Pt() / 1.e3 < 500.)
    updateFlag(m_vars.eventFlag, DiLeptonCuts::PtZLT500);
  else
    updateFlag(m_vars.eventFlag, DiLeptonCuts::PtZGT500);

  bool passPtvResolved(false), passPtvMerged(false);
  if (m_model == Model::AZh || m_model == Model::HVT) {
    checkPTVCut(m_vars.ZVec, m_vars.resolvedVH.vec_corr,
                m_vars.mergedVH.vec_corr, passPtvResolved, passPtvMerged);
    if (passPtvResolved)
      updateFlag(m_vars.eventFlag, DiLeptonCuts::PtZResolved);
    if (passPtvMerged) updateFlag(m_vars.eventFlag, DiLeptonCuts::PtZMerged);
  }

  //  **Selection** : C11-C14: jet multiplicity cuts
  //--------------------------------------------------
  if (m_physicsMeta.nJets == 2) {
    updateFlag(m_vars.eventFlag, DiLeptonCuts::Exactly2Jets);
  } else if (m_physicsMeta.nJets == 3) {
    updateFlag(m_vars.eventFlag, DiLeptonCuts::Exactly3Jets);
  } else if (m_physicsMeta.nJets == 4) {
    updateFlag(m_vars.eventFlag, DiLeptonCuts::Exactly4Jets);
  } else if (m_physicsMeta.nJets >= 5) {
    updateFlag(m_vars.eventFlag, DiLeptonCuts::AtLeast5Jets);
  }

  int nFatJetTags = 0;
  // Merged selection
  if (m_vars.fatJets.size() > 0) {
    nFatJetTags = Props::nBTags.get(m_vars.fatJets.at(0));
    // at least one track-jet in lead FJ
    updateFlag(m_vars.eventFlag, DiLeptonCuts::AtLeast1FJ);
    // count number of b-tagged leading track jets
    if (m_vars.trackJetsInLeadFJ.size() >= 1) {
      updateFlag(m_vars.eventFlag, DiLeptonCuts::LeadFJ1pTJ);
      if (nFatJetTags == 2)
        updateFlag(m_vars.eventFlag, DiLeptonCuts::LeadFJ2b);
    }
    if (m_vars.trackJetsInLeadFJ.size() > 1) {
      updateFlag(m_vars.eventFlag, DiLeptonCuts::LeadFJ2TJ);
    }
  }
  /////////////////////////////////////////////////////////////
  //  **Fill cutflow** :
  fill_2lepCutFlow(m_vars.eventFlag, m_physicsMeta.nJets, m_physicsMeta.nTags,
                   nFatJetTags, isMu, isE);

  ///////////////////////////////////////////////////
  // HISTO FILLING before categorisations        ///
  ///////////////////////////////////////////////////
  if (!m_doOnlyInputs) {
    EL_CHECK("AnalysisReader_VHreso2Lep::run_2Lep_analysis",
             fill_nJetHistos(m_vars.signalJets, "Sig"));
    if (!m_doReduceFillHistos)
      EL_CHECK("AnalysisReader_VHreso2Lep::run_2Lep_analysis",
               fill_nJetHistos(m_vars.forwardJets, "Fwd"));
    EL_CHECK("AnalysisReader_VHreso2Lep::run_2Lep_analysis",
             fill_nJetHistos(m_vars.fatJets, "Fat"));
    if (!m_doReduceFillHistos)
      EL_CHECK("AnalysisReader_VHreso2Lep::run_2Lep_analysis",
               fill_nJetHistos(m_vars.trackJets, "Trk"));
  }

  ///////////////////////////////////////////////////
  // DEFINITION OF CUTS FOR DIFFERENT REGIONS     ///
  ///////////////////////////////////////////////////

  // --- Cuts for resolved analysis --- //
  std::vector<unsigned long int> cuts_common_resolved = {
      DiLeptonCuts::Trigger,         DiLeptonCuts::LeadLepPt27,
      DiLeptonCuts::OSLeptons,       DiLeptonCuts::AtLeast2Jets,
      DiLeptonCuts::AtLeast2SigJets, DiLeptonCuts::PtB145};
  std::vector<unsigned long int> cuts_mBBcr_resolved_baseline = {
      DiLeptonCuts::SFLeptons, DiLeptonCuts::MllZwindowResolved,
      DiLeptonCuts::METHTResolved, DiLeptonCuts::PtZResolved};
  std::vector<unsigned long int> cuts_restrict_mBBcr_resolved = {
      DiLeptonCuts::mHCorrInsideSideBandResolved};
  std::vector<unsigned long int> cuts_SR_resolved = {
      DiLeptonCuts::mHCorrResolved};
  std::vector<unsigned long int> cuts_top_resolved;
  if (m_model == Model::AZh || m_model == Model::HVT)
    cuts_top_resolved = {DiLeptonCuts::MllZwindowResolved,
                         DiLeptonCuts::mHCorrResolved,
                         DiLeptonCuts::PtZResolved};
  else if (m_model == Model::CUT || m_model == Model::MVA)
    cuts_top_resolved = {DiLeptonCuts::MllAbove40};

  if (m_model == Model::CUT) {
    cuts_common_resolved = {DiLeptonCuts::LeadLepPt27,
                            DiLeptonCuts::Trigger,
                            DiLeptonCuts::OSLeptons,
                            DiLeptonCuts::AtLeast2Jets,
                            DiLeptonCuts::AtLeast2SigJets,
                            DiLeptonCuts::PtB145,
                            DiLeptonCuts::dRCut};
    cuts_mBBcr_resolved_baseline = {DiLeptonCuts::SFLeptons,
                                    DiLeptonCuts::MllZwindowResolved,
                                    DiLeptonCuts::METHTResolved};
    cuts_SR_resolved = {};
    cuts_top_resolved = {DiLeptonCuts::MllZwindowResolved};
  }

  if (m_model == Model::MVA) {
    cuts_common_resolved = {
        DiLeptonCuts::LeadLepPt27,     DiLeptonCuts::Trigger,
        DiLeptonCuts::OSLeptons,       DiLeptonCuts::AtLeast2Jets,
        DiLeptonCuts::AtLeast2SigJets, DiLeptonCuts::PtB145};
    cuts_mBBcr_resolved_baseline = {DiLeptonCuts::SFLeptons,
                                    DiLeptonCuts::MllZwindowResolved};
    cuts_SR_resolved = {};
    cuts_top_resolved = {DiLeptonCuts::MllZwindowResolved};
  }

  // mBBcr = common + mBBcr-baseline + restriction to 50-200 GeV
  std::vector<unsigned long int> cuts_mBBcr_resolved;
  cuts_mBBcr_resolved.insert(cuts_mBBcr_resolved.end(),
                             cuts_common_resolved.begin(),
                             cuts_common_resolved.end());
  cuts_mBBcr_resolved.insert(cuts_mBBcr_resolved.end(),
                             cuts_mBBcr_resolved_baseline.begin(),
                             cuts_mBBcr_resolved_baseline.end());
  cuts_mBBcr_resolved.insert(cuts_mBBcr_resolved.end(),
                             cuts_restrict_mBBcr_resolved.begin(),
                             cuts_restrict_mBBcr_resolved.end());
  // SR = (SR-specific) + common + mBBcr-baseline
  cuts_SR_resolved.insert(cuts_SR_resolved.end(), cuts_common_resolved.begin(),
                          cuts_common_resolved.end());
  cuts_SR_resolved.insert(cuts_SR_resolved.end(),
                          cuts_mBBcr_resolved_baseline.begin(),
                          cuts_mBBcr_resolved_baseline.end());
  // topCR = top-specific + common
  cuts_top_resolved.insert(cuts_top_resolved.end(),
                           cuts_common_resolved.begin(),
                           cuts_common_resolved.end());

  // --- Cuts for merged analysis --- //
  std::vector<unsigned long int> cuts_common_merged = {
      DiLeptonCuts::Trigger,         DiLeptonCuts::LeadLepPt27,
      DiLeptonCuts::AtLeast1FJ,      DiLeptonCuts::LeadFJ1pTJ,
      DiLeptonCuts::pTBothLeptons25, DiLeptonCuts::MuonEtaLT2p5};
  std::vector<unsigned long int> cuts_mBBcr_merged_baseline = {
      DiLeptonCuts::SFLeptons, DiLeptonCuts::MllZwindowMerged,
      DiLeptonCuts::METHTMerged, DiLeptonCuts::PtZMerged};
  std::vector<unsigned long int> cuts_restrict_mBBcr_merged = {
      DiLeptonCuts::mHCorrInsideSideBandMerged};
  std::vector<unsigned long int> cuts_SR_merged = {DiLeptonCuts::mHCorrMerged};
  std::vector<unsigned long int> cuts_top_merged;
  if (m_model == Model::AZh || m_model == Model::HVT)
    cuts_top_merged = {DiLeptonCuts::MllZwindowMerged,
                       DiLeptonCuts::mHCorrMerged, DiLeptonCuts::PtZMerged};
  else if (m_model == Model::CUT || m_model == Model::MVA)
    cuts_top_merged = {DiLeptonCuts::MllAbove40};

  if (m_model == Model::CUT) {
    cuts_common_merged = {DiLeptonCuts::LeadLepPt27, DiLeptonCuts::Trigger,
                          DiLeptonCuts::OSLeptons, DiLeptonCuts::AtLeast1FJ,
                          DiLeptonCuts::LeadFJ2TJ};
    cuts_mBBcr_merged_baseline = {DiLeptonCuts::SFLeptons,
                                  DiLeptonCuts::MllZwindowMerged,
                                  DiLeptonCuts::METHTMerged};
    cuts_SR_merged = {};
    cuts_top_merged = {DiLeptonCuts::MllZwindowMerged};
  }

  // mBBcr = common + mBBcr-baseline + restriction to 50-200 GeV
  std::vector<unsigned long int> cuts_mBBcr_merged;
  cuts_mBBcr_merged.insert(cuts_mBBcr_merged.end(), cuts_common_merged.begin(),
                           cuts_common_merged.end());
  cuts_mBBcr_merged.insert(cuts_mBBcr_merged.end(),
                           cuts_mBBcr_merged_baseline.begin(),
                           cuts_mBBcr_merged_baseline.end());
  cuts_mBBcr_merged.insert(cuts_mBBcr_merged.end(),
                           cuts_restrict_mBBcr_merged.begin(),
                           cuts_restrict_mBBcr_merged.end());
  // SR = (SR-specific) + common + mBBcr-baseline
  cuts_SR_merged.insert(cuts_SR_merged.end(), cuts_common_merged.begin(),
                        cuts_common_merged.end());
  cuts_SR_merged.insert(cuts_SR_merged.end(),
                        cuts_mBBcr_merged_baseline.begin(),
                        cuts_mBBcr_merged_baseline.end());
  // topCR = top-specific + common
  cuts_top_merged.insert(cuts_top_merged.end(), cuts_common_merged.begin(),
                         cuts_common_merged.end());

  /////////////////////////////////////////////////////////////
  //      DECISION TO USE MERGED OR RESOLVED ANALYSIS        //
  /////////////////////////////////////////////////////////////
  selectRegime(m_vars.eventFlag, m_vars.ZVec.Pt(), cuts_SR_resolved,
               cuts_SR_merged, cuts_mBBcr_resolved, cuts_mBBcr_merged,
               cuts_top_resolved, cuts_top_merged);

  /////////////////////////////////////////////////////////////
  //         EVENT CATEGORIZATION BASED ON BIT FLAG          //
  /////////////////////////////////////////////////////////////

  // choose right analysis model
  //---------------------------
  if (m_model == Model::MVA) {
    m_histNameSvc->set_analysisType(HistNameSvc::AnalysisType::MVA);
  } else if (m_model == Model::CUT) {
    m_histNameSvc->set_analysisType(HistNameSvc::AnalysisType::CUT);
  } else if (m_model == Model::AZh || m_model == Model::HVT) {
    m_histNameSvc->set_analysisType(HistNameSvc::AnalysisType::VHres);
  }

  // Event flavor, nTag, nJet, pTV
  //------------------------------
  if (m_physicsMeta.regime == PhysicsMetadata::Regime::resolved) {
    setevent_flavour(m_vars.selectedJets);
    m_histNameSvc->set_nTag(m_physicsMeta.nTags);
    m_histNameSvc->set_eventFlavour(m_physicsMeta.b1Flav, m_physicsMeta.b2Flav);
    if (!m_doMergePtVBins) m_histNameSvc->set_pTV(m_vars.ZVec.Pt());
    if (m_doMergeJetBins && m_physicsMeta.nJets >= 2)
      m_histNameSvc->set_nJet(-2);
    else
      m_histNameSvc->set_nJet(m_physicsMeta.nJets);

    if ((m_model == Model::CUT || m_model == Model::MVA)) {
      if (m_doMergeJetBins && m_physicsMeta.nJets >= 3) {
        m_histNameSvc->set_nJet(-3);
      } else if (m_physicsMeta.nJets >= 4) {
        m_histNameSvc->set_nJet(-4);
      } else {
        m_histNameSvc->set_nJet(m_physicsMeta.nJets);
      }
    }
  } else if (m_physicsMeta.regime == PhysicsMetadata::Regime::merged) {
    // Flavor labeling options - only 1 will be kept eventually
    if (m_truthLabeling == "TrackJetCone") {
      setevent_flavour(m_vars.trackjetsForLabelling);
    } else if (m_truthLabeling == "TrackJetGhostAssHadrons") {
      setevent_flavourGhost(m_vars.trackjetsForLabelling);
    } else if (m_truthLabeling == "FatJetGhostAssHadrons") {
      setevent_flavourGhost(m_vars.fatJets);
    } else if (m_truthLabeling == "TrackJetHybrid") {
      if (m_vars.trackJetsInLeadFJ.size() > 1)
        setevent_flavour(m_vars.trackjetsForLabelling);
      else
        setevent_flavourGhost(m_vars.trackjetsForLabelling);
    } else if (m_truthLabeling == "FatJetHybrid") {
      if (m_vars.trackJetsInLeadFJ.size() > 1)
        setevent_flavour(m_vars.trackjetsForLabelling);
      else
        setevent_flavourGhost(m_vars.fatJets);
    }

    if (m_doMergeModel)
      m_histNameSvc->set_nTrackJetInFatJet(m_vars.trackJetsInLeadFJ.size());
    m_histNameSvc->set_nTag(m_physicsMeta.nTagsInFJ);
    m_histNameSvc->set_nFatJet(m_vars.fatJets.size());
    if (!m_doMergePtVBins) m_histNameSvc->set_pTV(m_vars.ZVec.Pt());
    m_histNameSvc->set_nBTagTrackJetUnmatched(m_physicsMeta.nAddBTrkJets);
    m_histNameSvc->set_eventFlavour(m_physicsMeta.b1Flav, m_physicsMeta.b2Flav);
  }

  // Region
  //--------
  m_histNameSvc->set_description("");
  if (m_physicsMeta.regime == PhysicsMetadata::Regime::resolved) {
    if (passSpecificCuts(m_vars.eventFlag, cuts_SR_resolved)) {
      m_histNameSvc->set_description("SR");
    } else if (passSpecificCuts(m_vars.eventFlag, cuts_mBBcr_resolved) &&
               !passSpecificCuts(m_vars.eventFlag,
                                 {DiLeptonCuts::mHCorrResolved})) {
      if (!m_doMergeCR) {
        if (m_physicsMeta.mbbSideBandResolved ==
            PhysicsMetadata::MbbSideBandResolved::Low)
          m_histNameSvc->set_description("lowmBBcr");
        else if (m_physicsMeta.mbbSideBandResolved ==
                 PhysicsMetadata::MbbSideBandResolved::High)
          m_histNameSvc->set_description("highmBBcr");
      } else if (m_physicsMeta.mbbSideBandResolved ==
                     PhysicsMetadata::MbbSideBandResolved::Low ||
                 m_physicsMeta.mbbSideBandResolved ==
                     PhysicsMetadata::MbbSideBandResolved::High) {
        m_histNameSvc->set_description("mBBcr");
      }
    } else if (passSpecificCuts(m_vars.eventFlag, cuts_top_resolved) &&
               !passSpecificCuts(m_vars.eventFlag, {DiLeptonCuts::SFLeptons})) {
      m_histNameSvc->set_description("topemucr");
    }
  } else if (m_physicsMeta.regime == PhysicsMetadata::Regime::merged) {
    if (passSpecificCuts(m_vars.eventFlag, cuts_SR_merged)) {
      m_histNameSvc->set_description("SR");
    } else if (passSpecificCuts(m_vars.eventFlag, cuts_mBBcr_merged) &&
               !passSpecificCuts(m_vars.eventFlag,
                                 {DiLeptonCuts::mHCorrMerged})) {
      if (!m_doMergeCR) {
        if (m_physicsMeta.mbbSideBandMerged ==
            PhysicsMetadata::MbbSideBandMerged::Low)
          m_histNameSvc->set_description("lowmBBcr");
        else if (m_physicsMeta.mbbSideBandMerged ==
                 PhysicsMetadata::MbbSideBandMerged::High)
          m_histNameSvc->set_description("highmBBcr");
      } else if (m_physicsMeta.mbbSideBandMerged ==
                     PhysicsMetadata::MbbSideBandMerged::Low ||
                 m_physicsMeta.mbbSideBandMerged ==
                     PhysicsMetadata::MbbSideBandMerged::High) {
        m_histNameSvc->set_description("mBBcr");
      }
    } else if (passSpecificCuts(m_vars.eventFlag, cuts_top_merged) &&
               !passSpecificCuts(m_vars.eventFlag, {DiLeptonCuts::SFLeptons})) {
      m_histNameSvc->set_description("topemucr");
    }
  }

  /////////////////////////////////////////////////////////////
  // **RE-Definition** : Higgs, ZH - according to regime for HIST FILLING
  if (m_physicsMeta.regime == PhysicsMetadata::Regime::resolved) {
    m_vars.HVec = m_vars.resolvedH.vec_corr;
    m_vars.HVecPtReco = m_vars.resolvedH.vec_ptreco;
    if (m_physicsMeta.nJets > 2)
      m_vars.bbjVec = m_vars.resolvedH.vec_ptreco + m_vars.j3.vec;
    // if(m_histNameSvc->get_description() == "SR") m_vars.VHVec =
    // m_vars.resolvedVH.vec_resc;
    // else
    m_vars.VHVec = m_vars.resolvedVH.vec_corr;
  } else if (m_physicsMeta.regime == PhysicsMetadata::Regime::merged) {
    m_vars.HVec = m_vars.mergedH.vec_corr;
    if (m_vars.fatJets.size() > 1) {
      m_vars.bbjVec = m_vars.mergedH.vec_corr + m_vars.j2.vec_corr;
    }
    bool rescaleFatJetMass = false;
    m_config->getif<bool>("rescaleFatJetMass", rescaleFatJetMass);
    if (rescaleFatJetMass && m_histNameSvc->get_description() == "SR") {
      m_vars.VHVec = m_vars.mergedVH.vec_resc;
    } else {
      m_vars.VHVec = m_vars.mergedVH.vec_corr;
    }
  }

  if (m_physicsMeta.regime == PhysicsMetadata::Regime::resolved) {
    m_vars.HVec_gsc = m_vars.resolvedH.vec_gsc;
    m_vars.HVec_onemu = m_vars.resolvedH.vec_onemu;
    m_vars.HVec_ptreco = m_vars.resolvedH.vec_ptreco;
    m_vars.HVec_kf = m_vars.resolvedH.vec_kf;
    m_vars.VHVec_gsc = m_vars.resolvedVH.vec_gsc;
    m_vars.VHVec_onemu = m_vars.resolvedVH.vec_onemu;
    m_vars.VHVec_ptreco = m_vars.resolvedVH.vec_ptreco;
    m_vars.VHVec_kf = m_vars.resolvedVH.vec_kf;
    m_vars.lep1.vec_kf = m_vars.lep1.vec_kf_res;
    m_vars.lep2.vec_kf = m_vars.lep2.vec_kf_res;
  } else if (m_physicsMeta.regime == PhysicsMetadata::Regime::merged) {
    m_vars.HVec_gsc = m_vars.mergedH.vec_gsc;
    m_vars.HVec_onemu = m_vars.mergedH.vec_onemu;
    m_vars.HVec_ptreco = m_vars.mergedH.vec_ptreco;
    m_vars.HVec_kf = m_vars.mergedH.vec_kf;
    m_vars.VHVec_gsc = m_vars.mergedVH.vec_gsc;
    m_vars.VHVec_onemu = m_vars.mergedVH.vec_onemu;
    m_vars.VHVec_ptreco = m_vars.mergedVH.vec_ptreco;
    m_vars.VHVec_kf = m_vars.mergedVH.vec_kf;
    m_vars.lep1.vec_kf = m_vars.lep1.vec_kf_mer;
    m_vars.lep2.vec_kf = m_vars.lep2.vec_kf_mer;
  }

  /////////////////////////////////////////////////////////////
  // **Selection** : check if in blinded region
  // m_doBlinding = isBlindedRegion(m_vars.eventFlag, m_physicsMeta.regime ==
  // PhysicsMetadata::Regime::merged);
  m_doBlinding = false;

  /////////////////////////////////////////////////////////////
  // **Calculate** : b-tagging SF
  if (m_isMC) {
    if (m_physicsMeta.regime == PhysicsMetadata::Regime::resolved) {
      m_vars.btagWeight = computeBTagSFWeight(m_vars.signalJets,
                                              m_jetReader->getContainerName());
    } else if (m_physicsMeta.regime == PhysicsMetadata::Regime::merged) {
      m_vars.btagWeight = computeBTagSFWeight(
          m_vars.trackJetsForBTagging, m_trackJetReader->getContainerName());
    }
    m_weight *= m_vars.btagWeight;
    if (m_doTruthTagging)
      BTagProps::truthTagEventWeight.set(m_eventInfo, m_vars.btagWeight);
  }
  /////////////////////////////////////////////////////////////
  // store w/ / w/o PU weight as systematic variation
  // depends on what setting was chosen in config
  if (!m_config->get<bool>("applyPUWeight"))
    m_weightSysts.push_back({"_withPU", (float)(m_pileupReweight)});

  /////////////////////////////////////////////////////////////
  // Print event weights
  bool printEventWeightsForCutFlow = false;
  m_config->getif<bool>("printEventWeightsForCutFlow",
                        printEventWeightsForCutFlow);
  if (m_isMC && printEventWeightsForCutFlow &&
      passSpecificCuts(
          m_vars.eventFlag,
          {DiLeptonCuts::Trigger, DiLeptonCuts::SFLeptons,
           DiLeptonCuts::MllZwindowResolved, DiLeptonCuts::METHTResolved,
           DiLeptonCuts::AtLeast2Jets, DiLeptonCuts::AtLeast2SigJets,
           DiLeptonCuts::AtLeast1B, DiLeptonCuts::Exactly2B})) {
    printf(
        "Run: %7u  ,  Event: %8llu  -  PU = %6.3f  ,  Lepton = %5.3f  ,  "
        "Trigger = %5.3f , BTag = %5.3f\n",
        m_eventInfo->runNumber(), m_eventInfo->eventNumber(), m_pileupReweight,
        Props::leptonSF.get(m_eventInfo), m_vars.triggerSF, m_vars.btagWeight);
  }

  bool isNominal = false;
  if (m_currentVar == "Nominal") isNominal = true;
  /////////////////////////////////////////////////////////////
  // **CorrsAndSyst**:
  // resolved
  if (m_isMC && (m_model != Model::AZh && m_model != Model::HVT) &&
      (m_physicsMeta.regime == PhysicsMetadata::Regime::resolved) &&
      ((m_csCorrections.size() != 0) || (m_csVariations.size() != 0))) {
    // Check quantities (truth, reco, missing, ...) --> not well defined yet

    TLorentzVector j1Temp = m_vars.signalJets.at(0)->p4();
    TLorentzVector j2Temp = m_vars.signalJets.at(1)->p4();

    float cs_dr = fabs(j1Temp.DeltaR(j2Temp));
    float cs_dphi = fabs(j1Temp.DeltaPhi(j2Temp));
    float cs_vpt = m_vars.ZVec.Pt();
    float cs_mbb = (m_vars.j1.vec + m_vars.j2.vec)
                       .M();  // need to use uncorrected/scaled jets
    float cs_ptb1 = m_vars.j1.vec.Pt();
    float cs_ptb2 = m_vars.j2.vec.Pt();
    float cs_met = m_vars.met->met();
    float cs_njet = (m_physicsMeta.nSigJet <= 3 ? m_physicsMeta.nSigJet : 3);
    float cs_ntag = m_physicsMeta.nTags;
    float cs_truthPt = -1;
    float cs_avgTopPt = -1;

    // truth-level variables --> please check the definition if you want to
    // apply the systematic variations from CorrsAndSysts to samples different
    // from the nominal choices PowhegPythia6 for ttbar Sherpa for W+jets,
    // Z+jets, diboson
    truthVariablesCS(cs_truthPt, cs_avgTopPt);
    applyCS(cs_vpt, cs_mbb, cs_truthPt, cs_dphi, cs_dr, cs_ptb1, cs_ptb2,
            cs_met, cs_avgTopPt, cs_njet, cs_ntag);
  }

  // new Systematics for Vh resonance PRSR analysis
  bool nominalOnly = false;
  bool doModelSyst = true;
  m_config->getif<bool>("nominalOnly", nominalOnly);
  m_config->getif<bool>("doModelSyst", doModelSyst);

  if (m_isMC && (m_model == Model::AZh || m_model == Model::HVT)) {
    std::string Regime =
        (m_physicsMeta.regime == PhysicsMetadata::Regime::resolved) ? "Resolved"
                                                                    : "Boosted";
    bool isEMu =
        (passSpecificCuts(m_vars.eventFlag, cuts_top_resolved) &&
         !passSpecificCuts(m_vars.eventFlag, {DiLeptonCuts::SFLeptons}));

    if (m_histNameSvc->get_sample() == "ttbar") {
      if (!nominalOnly && isNominal && doModelSyst)
        fillTopSystslep_PRSR(Regime, m_physicsMeta.nAddBTrkJets,
                             m_vars.HVec.M(), m_vars.VHVec.M() / 1e3,
                             isEMu);  // mVH= GeV, mBB = MeV
    } else if (m_histNameSvc->get_sample() == "W" ||
               m_histNameSvc->get_sample() == "Wv22" ||
               m_histNameSvc->get_sample() == "Z" ||
               m_histNameSvc->get_sample() == "Zv22") {
      if (!nominalOnly && isNominal && doModelSyst)
        fillVjetsSystslep_PRSR(Regime, m_vars.HVec.M(),
                               m_vars.VHVec.M() / 1e3);  // mVH= GeV, mBB = MeV
    }
  }

  // Sherpa weight variations -> input for modelling syst
  // will only be run if appropriate syst are added to csVariations
  //(MUR0.5_MUF1_PDF261000 MUR2_MUF1_PDF261000 MUR1_MUF0.5_PDF261000
  // MUR1_MUF2_PDF261000 MUR0.5_MUF0.5_PDF261000 MUR2_MUF2_PDF261000
  // MUR1_MUF1_PDF25300 MUR1_MUF1_PDF13000)
  if (!nominalOnly && m_currentVar == "Nominal" && m_isMC)
    applySherpaVJet_EvtWeights();

  /////////////////////////////////////////////////////////////
  //          TREE and HISTO filling                         //
  /////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // **Trees**: only fill trees if there are more than 2 signal jets or more
  // than 1 fat-jet
  // easy tree
  EL_CHECK("run_2Lep_analysis()", FillEasyTree());

  // Fill the special tree for the orthogonality studies of the HVT combination
  EL_CHECK("run_2Lep_analysis()", fillOrthogonalityTree());
/*
  Long64_t m_EventNumber;
  m_EventNumber = m_eventInfo->eventNumber();
  m_MVAInputs2l->mBB = m_etree->get<float>("mBB");
  m_MVAInputs2l->dRBB = m_etree->get<float>("dRBBcor");
  m_MVAInputs2l->dEtaVBB = m_etree->get<float>("dEtaVBBPtReco");
  m_MVAInputs2l->dPhiVBB = m_etree->get<float>("dPhiVBBPtReco");
  m_MVAInputs2l->pTV = m_etree->get<float>("pTV");
  if (m_vars.j1.isTagged &&
      !m_vars.j2.isTagged) {  // 1 b-tagged => bjet is a first jet.
    m_MVAInputs2l->pTB1 = m_vars.j1.vec_ptreco.Pt();
    m_MVAInputs2l->pTB2 = m_vars.j2.vec_ptreco.Pt();
  } else if (!m_vars.j1.isTagged && m_vars.j2.isTagged) {
    m_MVAInputs2l->pTB1 = m_vars.j2.vec_ptreco.Pt();
    m_MVAInputs2l->pTB2 = m_vars.j1.vec_ptreco.Pt();
  } else {  // 0 or 2 b-tagged => order by their pT.
    if (m_vars.j1.vec.Pt() > m_vars.j2.vec.Pt() || !m_vars.j2.isTagged) {
      m_MVAInputs2l->pTB1 = m_vars.j1.vec_ptreco.Pt();
      m_MVAInputs2l->pTB2 = m_vars.j2.vec_ptreco.Pt();
    } else {
      m_MVAInputs2l->pTB1 = m_vars.j2.vec_ptreco.Pt();
      m_MVAInputs2l->pTB2 = m_vars.j1.vec_ptreco.Pt();
    }
  }
  m_MVAInputs2l->METSig =
      m_etree->get<float>("METHT") /
      sqrt(1000.);  // xml files assume that METHT is calculated in (GeV)
  m_MVAInputs2l->mLL = m_etree->get<float>("mLL");
  if (m_physicsMeta.nJets >= 3) {
    m_MVAInputs2l->mBBJ = m_vars.bbjVec.M();
    m_MVAInputs2l->pTJ3 = m_vars.j3.vec.Pt();
  }
  double mvaValue = -1;
  double mvaVZValue = -1;
  mvaValue *= 1.0;
  mvaVZValue *= 1.0;
  // VHbb
  if (m_physicsMeta.nJets == 2) {
    if (m_etree->get<float>("pTV") < 150.0 * 1000) {
      if (m_EventNumber % 2 == 0)
        mvaValue = m_MVAReaders2l2jLoB.at("mva2jLoPtB")->EvaluateMVA("BDT");
      else
        mvaValue = m_MVAReaders2l2jLoA.at("mva2jLoPtA")->EvaluateMVA("BDT");
    }
    if (m_etree->get<float>("pTV") >= 150.0 * 1000) {
      if (m_EventNumber % 2 == 0)
        mvaValue = m_MVAReaders2l2jHiB.at("mva2jHiPtB")->EvaluateMVA("BDT");
      else
        mvaValue = m_MVAReaders2l2jHiA.at("mva2jHiPtA")->EvaluateMVA("BDT");
    }
  } else if (m_physicsMeta.nJets >= 3) {
    if (m_etree->get<float>("pTV") < 150.0 * 1000) {
      if (m_EventNumber % 2 == 0)
        mvaValue = m_MVAReaders2l3jLoB.at("mva3jLoPtB")->EvaluateMVA("BDT");
      else
        mvaValue = m_MVAReaders2l3jLoA.at("mva3jLoPtA")->EvaluateMVA("BDT");
    }
    if (m_etree->get<float>("pTV") >= 150.0 * 1000) {
      if (m_EventNumber % 2 == 0)
        mvaValue = m_MVAReaders2l3jHiB.at("mva3jHiPtB")->EvaluateMVA("BDT");
      else
        mvaValue = m_MVAReaders2l3jHiA.at("mva3jHiPtA")->EvaluateMVA("BDT");
    }
  }
  // VZbb
  if (m_physicsMeta.nJets == 2) {
    if (m_etree->get<float>("pTV") < 150.0 * 1000) {
      if (m_EventNumber % 2 == 0)
        mvaVZValue =
            m_MVAReadersVZ2l2jLoB.at("mvaVZ2jLoPtB")->EvaluateMVA("BDT");
      else
        mvaVZValue =
            m_MVAReadersVZ2l2jLoA.at("mvaVZ2jLoPtA")->EvaluateMVA("BDT");
    }
    if (m_etree->get<float>("pTV") >= 150.0 * 1000) {
      if (m_EventNumber % 2 == 0)
        mvaVZValue =
            m_MVAReadersVZ2l2jHiB.at("mvaVZ2jHiPtB")->EvaluateMVA("BDT");
      else
        mvaVZValue =
            m_MVAReadersVZ2l2jHiA.at("mvaVZ2jHiPtA")->EvaluateMVA("BDT");
    }
  } else if (m_physicsMeta.nJets >= 3) {
    if (m_etree->get<float>("pTV") < 150.0 * 1000) {
      if (m_EventNumber % 2 == 0)
        mvaVZValue =
            m_MVAReadersVZ2l3jLoB.at("mvaVZ3jLoPtB")->EvaluateMVA("BDT");
      else
        mvaVZValue =
            m_MVAReadersVZ2l3jLoA.at("mvaVZ3jLoPtA")->EvaluateMVA("BDT");
    }
    if (m_etree->get<float>("pTV") >= 150.0 * 1000) {
      if (m_EventNumber % 2 == 0)
        mvaVZValue =
            m_MVAReadersVZ2l3jHiB.at("mvaVZ3jHiPtB")->EvaluateMVA("BDT");
      else
        mvaVZValue =
            m_MVAReadersVZ2l3jHiA.at("mvaVZ3jHiPtA")->EvaluateMVA("BDT");
    }
  }

  m_etree->SetBranchAndValue<float>("mva", mvaValue, -1);
  m_etree->SetBranchAndValue<float>("mvadiboson", mvaVZValue, -1);
*/
  double mvaValue = -999;
  double mvaVZValue = -999;
  // printMVA();  //print inputs and output of the MVA for debugging
  /////////////////////////////////////////////////////////////
  // **Histos**: only fill histos for defined regions
  bool GenerateSTXS = false;
  m_config->getif<bool>("GenerateSTXS", GenerateSTXS);
  if (!GenerateSTXS || m_physicsMeta.nTags > 2) {
    if (m_histNameSvc->get_description() != "") {
      if (!m_doOnlyInputs)
        EL_CHECK("AnalysisReader_VHreso2Lep::run_2Lep_analysis",
                 fill_jetHistos(m_vars.signalJets, m_vars.forwardJets));
      EL_CHECK("AnalysisReader_VHreso2Lep::run_2Lep_analysis",
               fill_jetSelectedHistos());
      if (!m_doOnlyInputs)
        EL_CHECK("AnalysisReader_VHreso2Lep::run_2Lep_analysis",
                 fill_fatJetHistos(m_vars.fatJets));
      EL_CHECK("AnalysisReader_VHreso2Lep::run_2Lep_analysis",
               fill_2LepHistos(isMu, isE));
    }
    if (m_histNameSvc->get_description() != "") {
      // m_histSvc->BookFillHist("mvadiboson", 1000, -1, 1, mvaVZValue, m_weight);
      // m_histSvc->BookFillHist("mva", 1000, -1, 1, mvaValue, m_weight);
      // m_histSvc->BookFillHist("mBBMVA", 500, 0., 500.,
      //                         m_MVAInputs2l->mBB / 1000., m_weight);
    }
  } else {
    if (STXS_GetBinFromEvtInfo() != -1 &&
        m_physicsMeta.nTags < 3) {  // check if it is the STXS signals
      STXS_FillYields();            // a 2-D histogram filled
      std::string temp_name = m_histNameSvc->get_sample();
      STXS_ReplaceName();  // replace the signal name with the STXS bin
      if (m_histNameSvc->get_description() != "") {
        if (!m_doOnlyInputs)
          EL_CHECK("AnalysisReader_VHreso2Lep::run_2Lep_analysis",
                   fill_jetHistos(m_vars.signalJets, m_vars.forwardJets));
        EL_CHECK("AnalysisReader_VHreso2Lep::run_2Lep_analysis",
                 fill_jetSelectedHistos());
        if (!m_doOnlyInputs)
          EL_CHECK("AnalysisReader_VHreso2Lep::run_2Lep_analysis",
                   fill_fatJetHistos(m_vars.fatJets));
        EL_CHECK("AnalysisReader_VHreso2Lep::run_2Lep_analysis",
                 fill_2LepHistos(isMu, isE));
      }
      if (m_histNameSvc->get_description() != "") {
        // m_histSvc->BookFillHist("mvadiboson", 1000, -1, 1, mvaVZValue,
        //                         m_weight);
        // m_histSvc->BookFillHist("mva", 1000, -1, 1, mvaValue, m_weight);
        // m_histSvc->BookFillHist("mBBMVA", 500, 0., 500.,
        //                         m_MVAInputs2l->mBB / 1000., m_weight);
      }
      // set the histogram name back from STXS addition
      m_histNameSvc->set_sample(temp_name);
    }
  }

  // Info("AnalysisReader_VHreso2Lep::run_2Lep_analysis()", "Region: %s, and
  // nTags: %d, njet: %d ", m_histNameSvc->get_description().c_str(),
  // m_physicsMeta.nTags, m_histNameSvc->get_nJet());

  return EL::StatusCode::SUCCESS;
}  // run_2Lep_analysis

EL::StatusCode AnalysisReader_VHreso2Lep::FillEasyTree() {
  // fill all EasyTree branches
  if (m_physicsMeta.nSigJet < 2 && m_vars.fatJets.size() < 1) {
    // fill only if we fullfil resolved or merged jet requirements
    return EL::StatusCode::SUCCESS;
  }
  // just used for shuffeling
  (*m_etree)["L1Vec"] = m_vars.lep1.vec;
  (*m_etree)["L2Vec"] = m_vars.lep2.vec;
  (*m_etree)["ZVec"] = m_vars.ZVec;
  (*m_etree)["HVec"] = m_vars.HVec;
  (*m_etree)["VHVec"] = m_vars.VHVec;
  if (m_physicsMeta.nJets > 2) (*m_etree)["bbjVec"] = m_vars.bbjVec;

  // written to output
  m_etree->SetBranchAndValue<std::string>("Description",
                                          m_histNameSvc->get_description(), "");
  m_etree->SetBranchAndValue<std::string>("Sample",
                                          m_histNameSvc->getFullSample(), "");
  m_etree->SetBranchAndValue<std::string>("EventFlavor",
                                          m_histNameSvc->getEventFlavour(), "");
  if (m_isMC)
    m_etree->SetBranchAndValue<float>("MCChannelNumber", m_mcChannel, -1);
  m_etree->SetBranchAndValue<float>("EventWeight", m_weight, -1);
  if (m_isMC)
    m_etree->SetBranchAndValue<float>(
        "MCEventWeight", Props::MCEventWeight.get(m_eventInfo), -1);
  m_etree->SetBranchAndValue<float>("LumiWeight",
                                    Props::LumiWeight.get(m_eventInfo), -1);
  float NTruthJetWeight = Props::NTruthWZJets20.exists(m_eventInfo)
                              ? Props::NTruthWZJets20.get(m_eventInfo)
                              : -1;
  m_etree->SetBranchAndValue<float>("NTruthJetWeight", NTruthJetWeight, -1);
  m_etree->SetBranchAndValue<float>("PUWeight", m_pileupReweight, -1);
  m_etree->SetBranchAndValue<float>("BTagSF", m_vars.btagWeight, -1);
  m_etree->SetBranchAndValue<float>("TriggerSF", m_vars.triggerSF, -1);
  m_etree->SetBranchAndValue<float>("LeptonSF",
                                    Props::leptonSF.get(m_eventInfo), -1);
  m_etree->SetBranchAndValue<float>("JVTWeight", m_vars.JVTweight, -1);
  m_etree->SetBranchAndValue<int>("EventNumber", m_eventInfo->eventNumber(),
                                  -1);
  m_etree->SetBranchAndValue<int>(
      "passedTrigger",
      passSpecificCuts(m_vars.eventFlag, {DiLeptonCuts::Trigger}), -1);
  m_etree->SetBranchAndValue<int>(
      "passedTauVeto",
      passSpecificCuts(m_vars.eventFlag, {DiLeptonCuts::tauVeto}), -1);
  m_etree->SetBranchAndValue<float>("AverageMu", m_averageMu, -1);
  m_etree->SetBranchAndValue<float>("ActualMu", m_actualMu, -1);
  m_etree->SetBranchAndValue<float>("AverageMuScaled", m_averageMuScaled, -1);
  m_etree->SetBranchAndValue<float>("ActualMuScaled", m_actualMuScaled, -1);
  int Nvtx = 0;
  if (Props::NVtx2Trks.exists(m_eventInfo)) {
    Nvtx = Props::NVtx2Trks.get(m_eventInfo);
  } else if (Props::NVtx3Trks.exists(m_eventInfo)) {
    Nvtx = Props::NVtx3Trks.get(m_eventInfo);
  }
  m_etree->SetBranchAndValue<int>("Nvtx", Nvtx, -1);
  m_etree->SetBranchAndValue<float>("ZPV", Props::ZPV.get(m_eventInfo), -1);
  m_etree->SetBranchAndValue<float>("ptL1", m_vars.lep1.vec.Pt(), -1);
  m_etree->SetBranchAndValue<float>("ptL2", m_vars.lep2.vec.Pt(), -1);
  m_etree->SetBranchAndValue<float>("etaL1", m_vars.lep1.vec.Eta(), -1);
  m_etree->SetBranchAndValue<float>("etaL2", m_vars.lep2.vec.Eta(), -1);
  m_etree->SetBranchAndValue<float>("flavL1", m_vars.lep1.flav, -1);
  m_etree->SetBranchAndValue<float>("flavL2", m_vars.lep2.flav, -1);
  int isMediumL1 = 1;
  int isMediumL2 = 1;
  if (m_vars.el1 != nullptr) {
    isMediumL1 = Props::isMediumLH.get(m_vars.el1);
  }
  if (m_vars.el2 != nullptr) {
    isMediumL2 = Props::isMediumLH.get(m_vars.el2);
  }


/*
  if (m_vars.mu1)
  {
    cout << m_vars.mu1->auxdata<float>("d0sigBL") << endl;
    m_etree->SetBranchAndValue<float>("d0sigBL1", m_vars.lep1.vec.Phi(), -1);
  }
  if (m_vars.mu2)
  {
    cout << m_vars.mu2->auxdata<float>("d0sigBL2") << endl;
  }
  if (m_vars.el1)
  {
    cout << m_vars.el1->auxdata<float>("d0sigBL1") << endl;
  }
  if (m_vars.el2)
  {
    cout << m_vars.el2->auxdata<float>("d0sigBL2") << endl;
  }
*/

  m_etree->SetBranchAndValue<float>("d0sigBLL1", m_vars.lep1.d0sigBL, -1);
  m_etree->SetBranchAndValue<float>("d0sigBLL2", m_vars.lep2.d0sigBL, -1);
  m_etree->SetBranchAndValue<float>("z0sinThetaL1", m_vars.lep1.z0sinTheta, -1);
  m_etree->SetBranchAndValue<float>("z0sinThetaL2", m_vars.lep2.z0sinTheta, -1);

  m_etree->SetBranchAndValue<float>("PhiL1", m_vars.lep1.vec.Phi(), -1);
  m_etree->SetBranchAndValue<float>("PhiL2", m_vars.lep2.vec.Phi(), -1);
  //m_etree->SetBranchAndValue<float>("d0L2", m_vars.lep2.d0wrtPrimVtx(), -1);

  //m_etree->SetBranchAndValue<float>("d0sigBL", m_vars.lep1.d0sigBL(), -1);


  m_etree->SetBranchAndValue<int>("isMediumL1", isMediumL1, -1);
  m_etree->SetBranchAndValue<int>("isMediumL2", isMediumL2, -1);
  m_etree->SetBranchAndValue<float>("chargeL1", m_vars.lep1.charge, -1);
  m_etree->SetBranchAndValue<float>("chargeL2", m_vars.lep2.charge, -1);
  m_etree->SetBranchAndValue<float>(
      "dEtaLL", fabs(m_vars.lep1.vec.Eta() - m_vars.lep2.vec.Eta()), -1);
  m_etree->SetBranchAndValue<float>(
      "dPhiLL", fabs(m_vars.lep1.vec.DeltaPhi(m_vars.lep2.vec)), -1);
  m_etree->SetBranchAndValue<float>("mLL", m_vars.ZVec.M(), -1);
  m_etree->SetBranchAndValue<float>("pTV", m_vars.ZVec.Pt(), -1);
  m_etree->SetBranchAndValue<float>("MET", m_vars.met->met(), -1);
  m_etree->SetBranchAndValue<float>("HT", m_vars.HTsmallR, -1);
  m_etree->SetBranchAndValue<float>("METHT", m_vars.METHTsmallR, -1);
  m_etree->SetBranchAndValue<float>(
      "dEtaBBcor", fabs(m_vars.j1.vec_corr.Eta() - m_vars.j2.vec_corr.Eta()),
      -1);
  m_etree->SetBranchAndValue<float>(
      "dPhiBBcor", fabs(m_vars.j1.vec_corr.DeltaPhi(m_vars.j2.vec_corr)), -1);
  m_etree->SetBranchAndValue<float>(
      "dRBBcor", fabs(m_vars.j1.vec_corr.DeltaR(m_vars.j2.vec_corr)), -1);
  m_etree->SetBranchAndValue<float>(
      "dPhiVBB", fabs(m_vars.ZVec.DeltaPhi(m_vars.HVec)), -1);
  m_etree->SetBranchAndValue<float>(
      "dEtaVBB", fabs(m_vars.ZVec.Eta() - m_vars.HVec.Eta()), -1);
  m_etree->SetBranchAndValue<float>(
      "dPhiVBBPtReco", fabs(m_vars.ZVec.DeltaPhi(m_vars.HVecPtReco)), -1);
  m_etree->SetBranchAndValue<float>(
      "dEtaVBBPtReco", fabs(m_vars.ZVec.Eta() - m_vars.HVecPtReco.Eta()), -1);
  m_etree->SetBranchAndValue<float>("dRVBB",
                                    fabs(m_vars.ZVec.DeltaR(m_vars.HVec)), -1);
  m_etree->SetBranchAndValue<float>("pTVH", m_vars.VHVec.Pt(), -1);
  m_etree->SetBranchAndValue<float>("mVH", m_vars.VHVec.M(), -1);
  m_etree->SetBranchAndValue<float>("mVHres", m_vars.resolvedVH.vec_resc.M(),
                                    -1);
  m_etree->SetBranchAndValue<float>("mVHmerg", m_vars.mergedVH.vec_corr.M(),
                                    -1);
  m_etree->SetBranchAndValue<float>("mBBres", m_vars.resolvedH.vec_corr.M(),
                                    -1);
  m_etree->SetBranchAndValue<float>("mBBmerg", m_vars.mergedH.vec_corr.M(), -1);

  TLorentzVector mLB =
      defineMinDRLepBSyst(m_vars.lep1, m_vars.lep2, m_vars.j1, m_vars.j2);
  m_etree->SetBranchAndValue<float>("mLB", mLB.M(), -1);
  m_etree->SetBranchAndValue<int>("nTaus", m_vars.taus.size(), -1);
  // m_etree->SetBranchAndValue("blind", (), false);
  EL_CHECK("AnalysisReader_VHreso2Lep::run_2Lep_analysis",
           compute_jetSelectedQuantities(m_vars.selectedJets));
  // Adding here the fat-jet and track-jet related quantities in order to
  // avoid changing compute_jetSelectedQuantities which is used by the other
  // analyses
  m_etree->SetBranchAndValue<int>("nFatJets", m_vars.fatJets.size(), -1);
  if (m_vars.fatJets.size() >= 1) {
    m_etree->SetBranchAndValue<int>("nbTagsInFJ", m_physicsMeta.nTagsInFJ, -1);
    m_etree->SetBranchAndValue<int>(
        "nTrkjetsInFJ", Props::nTrackJets.get(m_vars.fatJets.at(0)), -1);
  }
  m_etree->SetBranchAndValue<int>("nbTagsOutsideFJ", m_physicsMeta.nAddBTrkJets,
                                  -1);
  m_etree->SetBranchAndValue<int>("nTrkJets", m_vars.trackJets.size(), -1);
  m_etree->SetBranchAndValue<int>("nBTrkJets", m_vars.nTaggedTrkJets, -1);
  m_etree->SetBranchAndValue<float>("HTfat", m_vars.HTlargeR, -1);
  m_etree->SetBranchAndValue<float>("METHTfat", m_vars.METHTlargeR, -1);
  if (m_vars.fatJets.size() >= 1)
    (*m_etree)["fatj1Vec"] = m_vars.fatJets.at(0)->p4();
  if (m_vars.fatJets.size() >= 2)
    (*m_etree)["fatj2Vec"] = m_vars.fatJets.at(1)->p4();
  if (m_vars.fatJets.size() >= 3)
    (*m_etree)["fatj3Vec"] = m_vars.fatJets.at(2)->p4();
  if (m_vars.trackJets.size() >= 1)
    (*m_etree)["trkj1Vec"] = m_vars.trackJets.at(0)->p4();
  if (m_vars.trackJets.size() >= 2)
    (*m_etree)["trkj2Vec"] = m_vars.trackJets.at(1)->p4();
  if (m_vars.trackJets.size() >= 3)
    (*m_etree)["trkj3Vec"] = m_vars.trackJets.at(2)->p4();
  if (m_vars.trackJetsInLeadFJ.size() >= 1) {
    (*m_etree)["trkj1LeadFJVec"] = m_vars.trackJetsInLeadFJ.at(0)->p4();
    // m_etree->SetBranchAndValue<float>("trkj1LeadFJ_MV2c20",
    // Props::MV2c20.get(trackJetsInLeadFJ.at(0)), -1);
    if (m_vars.trackJetsInLeadFJ.size() >= 2) {
      (*m_etree)["trkj2LeadFJVec"] = m_vars.trackJetsInLeadFJ.at(1)->p4();
      // m_etree->SetBranchAndValue<float>("trkj2LeadFJ_MV2c20",
      // Props::MV2c20.get(trackJetsInLeadFJ.at(1)), -1);
    }
  }

  m_etree->SetBranchAndValue<float>("GSCMbb", m_vars.HVec_gsc.M(), -1);
  m_etree->SetBranchAndValue<float>("OneMuMbb", m_vars.HVec_onemu.M(), -1);
  m_etree->SetBranchAndValue<float>("PtRecoMbb", m_vars.HVec_ptreco.M(), -1);
  m_etree->SetBranchAndValue<float>("KFMbb", m_vars.HVec_kf.M(), -1);

  m_etree->SetBranchAndValue<float>("GSCMvh", m_vars.VHVec_gsc.M(), -1);
  m_etree->SetBranchAndValue<float>("OneMuMvh", m_vars.VHVec_onemu.M(), -1);
  m_etree->SetBranchAndValue<float>("PtRecoMvh", m_vars.VHVec_ptreco.M(), -1);
  m_etree->SetBranchAndValue<float>("KFMvh", m_vars.VHVec_kf.M(), -1);

  // Regime
  string Regime = "none";
  if (m_physicsMeta.regime == PhysicsMetadata::Regime::merged) {
    Regime = "merged";
  } else if (m_physicsMeta.regime == PhysicsMetadata::Regime::resolved)
    Regime = "resolved";

  m_etree->SetBranchAndValue<std::string>("Regime", Regime, "none");

  if (m_writeEasyTree && m_histNameSvc->get_isNominal()) {
    m_etree->Fill();
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode AnalysisReader_VHreso2Lep::setLeptonVariables() {
  // overloads base clas method using member variables
  AnalysisReader_VHQQ2Lep::setLeptonVariables(
      m_vars.lep1, m_vars.lep2, m_vars.el1, m_vars.el2, m_vars.mu1, m_vars.mu2);

  return EL::StatusCode::SUCCESS;
}
EL::StatusCode AnalysisReader_VHreso2Lep::setTwoLeptonFlavour() {
  // overloads base clas method using member variables
  AnalysisReader_VHQQ2Lep::setTwoLeptonFlavour(m_vars.lep1, m_vars.lep2);

  return EL::StatusCode::SUCCESS;
}
EL::StatusCode AnalysisReader_VHreso2Lep::fillOrthogonalityTree() {
  double mass;
  std::vector<int> regions;
  // The convention for the regions is defined in this Google Doc
  // https://docs.google.com/spreadsheets/d/1_N9L0z5L2DNGxt2NXsAH-O5Owykt7S7dKFNX5Y2h2To/edit#gid=0
  // We have the following categories
  //   - 551: resolved 1 tag
  //   - 552: resolved 2 tag
  //   - 561: merged 1 tag
  //   - 562: merged 2 tag
  // By convention a vector of regions that the event passes is written to the
  // EasyTree. VHreso SRs are disjoint, so the vector will have only a single
  // entry.
  // if (m_analysisStrategy != "Resolved" && m_analysisStrategy != "") {
  //   Error("fillOrthogonalityTree",
  //         "Not implemented for analysisStrategy != Resolved");
  //   return EL::StatusCode::FAILURE;
  // }

  // FIXME RLes Copied from above, should probably store as a flag
  std::vector<unsigned long int> cuts_common_merged = {
      DiLeptonCuts::Trigger,         DiLeptonCuts::LeadLepPt27,
      DiLeptonCuts::AtLeast1FJ,      DiLeptonCuts::LeadFJ1pTJ,
      DiLeptonCuts::pTBothLeptons25, DiLeptonCuts::MuonEtaLT2p5};
  std::vector<unsigned long int> cuts_mBBcr_merged_baseline = {
      DiLeptonCuts::SFLeptons, DiLeptonCuts::MllZwindowMerged,
      DiLeptonCuts::METHTMerged, DiLeptonCuts::PtZMerged};
  std::vector<unsigned long int> cuts_common_resolved = {
      DiLeptonCuts::Trigger,         DiLeptonCuts::LeadLepPt27,
      DiLeptonCuts::OSLeptons,       DiLeptonCuts::AtLeast2Jets,
      DiLeptonCuts::AtLeast2SigJets, DiLeptonCuts::PtB145};
  std::vector<unsigned long int> cuts_mBBcr_resolved_baseline = {
      DiLeptonCuts::SFLeptons, DiLeptonCuts::MllZwindowResolved,
      DiLeptonCuts::METHTResolved, DiLeptonCuts::PtZResolved};

  std::vector<unsigned long int> cuts_SR_merged = {DiLeptonCuts::mHCorrMerged};
  cuts_SR_merged.insert(cuts_SR_merged.end(), cuts_common_merged.begin(),
                        cuts_common_merged.end());
  cuts_SR_merged.insert(cuts_SR_merged.end(),
                        cuts_mBBcr_merged_baseline.begin(),
                        cuts_mBBcr_merged_baseline.end());

  std::vector<unsigned long int> cuts_SR_resolved = {
      DiLeptonCuts::mHCorrResolved};
  cuts_SR_resolved.insert(cuts_SR_resolved.end(), cuts_common_resolved.begin(),
                          cuts_common_resolved.end());
  cuts_SR_resolved.insert(cuts_SR_resolved.end(),
                          cuts_mBBcr_resolved_baseline.begin(),
                          cuts_mBBcr_resolved_baseline.end());

  bool passesResolvedSR = passSpecificCuts(
      m_vars.eventFlag,
      cuts_SR_resolved);  // FIXME RLES is this the proper defintion of SR?
  bool passesMergedSR = passSpecificCuts(m_vars.eventFlag, cuts_SR_merged);
  int region = 0;
  int ntags = 0;
  if (passesResolvedSR) {
    region = 550;
    ntags = m_physicsMeta.nTags;
    mass = m_vars.resolvedVH.vec_corr.M();
  } else if (passesMergedSR) {
    region = 560;
    ntags = m_physicsMeta.nTagsInFJ;
    mass = m_vars.mergedVH.vec_corr.M();
  } else {
    if (m_debug) {
      Info("fillOrthogonalityTree()", "Event is not in a SR. Returning.");
    }
    return EL::StatusCode::SUCCESS;
  }
  if (ntags != 1 && ntags != 2) {
    if (m_debug) {
      Info("fillOrthogonalityTree()",
           "Event has neither 1 nor 2 tags. Returning.");
    }
    return EL::StatusCode::SUCCESS;
  }
  AnalysisReader_VHQQ::fillOrthogonalityTree(mass, "mVH", {region + ntags});
  return EL::StatusCode::SUCCESS;
}
