// -*- C++ -*-
//
// Package:    Ho_MuGunAnalyzer
// Class:      Ho_MuGunAnalyzer
// 
/**\class Ho_MuGunAnalyzer Ho_MuGunAnalyzer.cc ho_mugun/Ho_MuGunAnalyzer/src/Ho_MuGunAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yusuf Erdogan,32 4-B20,+41227676487
//         Created:  Fri Aug 16 09:54:59 CEST 2013
// $Id$
//
//


// system include files
#include <stdio.h> 
#include <memory>
#include <math.h> 

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/RefToBaseVector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"


#include "RecoMuon/MuonIdentification/interface/MuonHOAcceptance.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputer.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputerRcd.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrack.h"

#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetMatchInfo.h"

#include "TH2F.h"
#include "TGraph.h"
#include "TMultiGraph.h"
//
// class declaration
//

class Ho_MuGunAnalyzer : public edm::EDAnalyzer {
   public:
      explicit Ho_MuGunAnalyzer(const edm::ParameterSet&);
      ~Ho_MuGunAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      
      std::map<std::string, TH1*> histos_;
	  edm::Service<TFileService> theFileService;
	  
	  TFile* theOutputFile;
	  TGraph* etaphiScatter_ge;
	  TGraph* etaphiScatter_be;
	  TGraph* etaphiScatter_ne;
	  
	  TMultiGraph * HODeadRegions;
	  TMultiGraph * HOSiPMRegions;
      
      MuonHOAcceptance* theMuonHOAcceptance;
      //TrackAssociatorParameters assocParams;
      //TrackDetectorAssociator assoc;
      const CaloGeometry *theGeometry;
      
      int i,j,k;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Ho_MuGunAnalyzer::Ho_MuGunAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   
   	//assoc.useDefaultPropagator();
   	
   	//assocParams = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
   
   i = 0;
   j = 0;
   k = 0;
   
	histos_["gen_p"] = theFileService->make<TH1F> ("gen_p", "gen_p", 100, 0, 200.);
	histos_["gen_phi"] = theFileService->make<TH1F> ("gen_phi", "gen_phi", 80, -4., 4.);
	histos_["gen_eta"] = theFileService->make<TH1F> ("gen_eta", "gen_eta", 30, -1.5, 1.5);
	histos_["gen_pt"] = theFileService->make<TH1F> ("gen_pt", "gen_pt", 100, 0, 200.);
	
	histos_["gen_id"] = theFileService->make<TH1F> ("gen_id", "gen_id", 200, -100, 100);
	
	histos_["hit_freq"] = theFileService->make<TH1F> ("hit_freq", "hit_freq", 301, -0.5, 300.5);
	
	//histos_["hit_ids"]
	
	histos_["part_geneta_corr"] = theFileService->make<TH1F> ("part_geneta_corr", "part_geneta_corr", 300, -1.5, 1.5);
	histos_["whole_geneta_nocorr"] = theFileService->make<TH1F> ("whole_geneta_nocorr", "whole_geneta_nocorr", 300, -1.5, 1.5);
	histos_["part_geneta_nocorr"] = theFileService->make<TH1F> ("part_geneta_nocorr", "part_geneta_nocorr", 300, -1.5, 1.5);
	
	histos_["geneta_aftercorr"] = theFileService->make<TH1F> ("geneta_aftercorr", "geneta_aftercorr", 400, -0.4, 0.4);
	
	histos_["ieta_iphi_hit"] = theFileService->make<TH2F> ("ieta_iphi_hit", "ieta_iphi_hit", 32, -16, 16, 80, 0, 80);
	
	histos_["hit_energy"] = theFileService->make<TH1F> ("hit_energy", "hit_energy", 100, 0., 10.);
	
	histos_["hit_maxenergy_wheel0"] = theFileService->make<TH1F> ("hit_maxenergy_wheel0", "hit_maxenergy_wheel0", 100, 0., 10.);
	histos_["hit_maxenergy_outer"] = theFileService->make<TH1F> ("hit_maxenergy_outer", "hit_maxenergy_outer", 100, 0., 10.);
	
	histos_["sumenergy"] = theFileService->make<TH1F> ("sumenergy", "sumenergy", 100, 0., 10.);
	histos_["sumenergy_wheel0"] = theFileService->make<TH1F> ("sumenergy_wheel0", "sumenergy_wheel0", 100, 0., 10.);
	histos_["sumenergy_outer"] = theFileService->make<TH1F> ("sumenergy_outer", "sumenergy_outer", 100, 0., 10.);
	
	histos_["gen_freq"] = theFileService->make<TH1F> ("gen_freq", "gen_freq", 2, -0.5, 1.5);
	histos_["gens_ingeomAccept"] = theFileService->make<TH1F> ("gens_ingeomAccept", "gens_ingeomAccept", 2, -0.5, 1.5);
	histos_["gens_inNotDead"] = theFileService->make<TH1F> ("gens_inNotDead", "gens_inNotDead", 2, -0.5, 1.5);
	histos_["gens_inSipm"] = theFileService->make<TH1F> ("gens_inSipm", "gens_inSipm", 2, -0.5, 1.5);
	
	histos_["hitvalidity"] = theFileService->make<TH1F> ("hitvalidity", "hitvalidity", 2, -0.5, 1.5);
	
	histos_["a1"] = theFileService->make<TH1F> ("a1", "a1", 2, -0.5, 1.5);
	histos_["a2"] = theFileService->make<TH1F> ("a2", "a2", 2, -0.5, 1.5);
	histos_["a3"] = theFileService->make<TH1F> ("a3", "a3", 2, -0.5, 1.5);
	histos_["a4"] = theFileService->make<TH1F> ("a4", "a4", 2, -0.5, 1.5);
	histos_["a5"] = theFileService->make<TH1F> ("a5", "a5", 2, -0.5, 1.5);
	histos_["a6"] = theFileService->make<TH1F> ("a6", "a6", 2, -0.5, 1.5);
	
	HODeadRegions = new TMultiGraph();
	HOSiPMRegions = new TMultiGraph();
	
	etaphiScatter_ge = new TGraph();
	etaphiScatter_be = new TGraph();
	etaphiScatter_ne = new TGraph();
	
	theOutputFile=new TFile("myFile.root","RECREATE");
	theOutputFile->Append(HODeadRegions);
	theOutputFile->Append(HOSiPMRegions);
}


Ho_MuGunAnalyzer::~Ho_MuGunAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Ho_MuGunAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

	//ESHandle < CaloGeometry > caloGeom;
	//iSetup.get<CaloGeometryRecord> ().get(caloGeom);

	//theGeometry = caloGeom.product();
	
	if (!MuonHOAcceptance::Inited()) MuonHOAcceptance::initIds(iSetup);
	
    Handle<std::vector<reco::GenParticle>> genparts;

  	iEvent.getByLabel("genParticles", genparts);
  	
  	reco::GenParticleCollection::const_iterator genPart = genparts->begin();

  	float genp = genPart->p();
  	float geneta = genPart->eta();
  	float genphi = genPart->phi();
  	float genpt = genPart->pt();
  	int genid = genPart->pdgId();
  	
  	//if(!(fabs(geneta) < 0.4)) return;
  	
  	histos_["gen_p"]->Fill(genp);
	histos_["gen_eta"]->Fill(geneta);
	histos_["gen_phi"]->Fill(genphi);
	histos_["gen_pt"]->Fill(genpt);
	histos_["gen_id"]->Fill(genid);
	
	histos_["gen_freq"]->Fill(genparts->size());
	
	bool inAccept = 1;
	bool inNotDead = 1;
	bool inSipm = 1;
	
	HODeadRegions = theMuonHOAcceptance->graphDeadRegions();
	HOSiPMRegions = theMuonHOAcceptance->graphSiPMRegions();
	
	// analyze only muons in geom.accept of HO
	if (!theMuonHOAcceptance->inGeomAccept(geneta, genphi, 0.04, 0.017)) inAccept = 0;
	//if (!theMuonHOAcceptance->inGeomAccept(geneta, genphi, 0., 0.)) inAccept = 0;
	// analyze only muons in non dead cells of HO
	if (!theMuonHOAcceptance->inNotDeadGeom(geneta, genphi, 0.04, 0.017)) inNotDead = 0;
	//if (!theMuonHOAcceptance->inNotDeadGeom(geneta, genphi, 0., 0.)) inNotDead = 0;
	// analyze only muons in SiPM tiles of HO
	//if (!theMuonHOAcceptance->inSiPMGeom(geneta, genphi, 0., 0.)) inSipm = 0;

	histos_["gens_ingeomAccept"]->Fill(inAccept);
	histos_["gens_inNotDead"]->Fill(inNotDead);
	histos_["gens_inSipm"]->Fill(inSipm);

	bool a1 = 0;
	bool a2 = 0;
	bool a3 = 0;
	bool a4 = 0;
	bool a5 = 0;
	bool a6 = 0;

	if(inAccept && !inNotDead && !inSipm) a1 = 1;
	histos_["a1"]->Fill(a1);
	if(inAccept && inNotDead && inSipm) a2 = 1;
	histos_["a2"]->Fill(a2);
	if(!inAccept && inNotDead && !inSipm) a3 = 1;
	histos_["a3"]->Fill(a3);
	if(!inAccept && inNotDead && inSipm) a4 = 1;
	histos_["a4"]->Fill(a4);
	if(!inAccept && !inNotDead && inSipm) a5 = 1;
	histos_["a5"]->Fill(a5);
	if(inAccept && !inNotDead && inSipm) a6 = 1;
	histos_["a6"]->Fill(a6);
	
	if (!(inAccept == 1 && inNotDead == 1 && inSipm == 1)) return;

	//GlobalPoint genvertex(genPart->vx(), genPart->vy(), genPart->vz());
 	//GlobalVector genmomentum(genPart->px(), genPart->py(), genPart->pz());
 	//TrackDetMatchInfo genMatch = assoc.associate(iEvent, iSetup, genmomentum, genvertex, genPart->charge(), assocParams);

	Handle<PCaloHitContainer> HcalHits;
	iEvent.getByLabel("g4SimHits", "HcalHits", HcalHits);
	
	histos_["hitvalidity"]->Fill(HcalHits.isValid());
	
	if (!HcalHits.isValid()) return;
	
	histos_["hit_freq"]->Fill(HcalHits->size());
	
	PCaloHitContainer::const_iterator simhit;
	
	std::map<uint32_t,double> homap;
	std::vector<double> hoTimeVec;
	std::map<uint32_t,double>::iterator hitsum;
	
	homap.clear();
	hoTimeVec.clear();
	
	for (simhit = HcalHits->begin(); simhit != HcalHits->end(); ++simhit) {
    	HcalDetId det(simhit->id());
    	
    	if (!(det.subdet() == 3)) continue;

    	std::map<uint32_t,double>::iterator lb = homap.lower_bound(simhit->id());
    	if ((lb != homap.end()) && (simhit->id() == lb->first)) {
      		homap[simhit->id()] += simhit->energy();
      	} else {
      		lb = homap.insert(lb, std::map<uint32_t,double>::value_type(simhit->id(),simhit->energy()));
    	}
	}

	//angle for the free-path-length correction
	double gentheta = -99.;
	double angle_rad = -99.;
	gentheta = 2*atan(exp(-1*geneta));
	angle_rad = 1.5708 - gentheta;
	
  	double sumsimhits = 0.;
  	double hit_maxenergy = -99.;
  	double maxEta = -99.;
  	double maxPhi = -99.;
  	
	for (hitsum = homap.begin(); hitsum != homap.end(); ++hitsum) {
		HcalDetId tmpId(hitsum->first);
      	sumsimhits += hitsum->second;
      	
      	if ((tmpId.ieta() == maxEta) && (tmpId.iphi() == maxPhi)) {
			hit_maxenergy += hitsum->second;
      	}
      	
      	if (hitsum->second > hit_maxenergy) {
			hit_maxenergy = hitsum->second;
			maxEta = tmpId.ieta();
			maxPhi = tmpId.iphi();
      	}
	}

  	histos_["sumenergy"]->Fill(sumsimhits*1000*cos(angle_rad));
	
	//simhit energies GeV -> MeV
	hit_maxenergy = hit_maxenergy*1000;
	
	double hit_maxenergy_corr = -999999;
	hit_maxenergy_corr = hit_maxenergy*cos(angle_rad);
	
	if(fabs(geneta) < 0.307) {
		histos_["hit_maxenergy_wheel0"]->Fill(hit_maxenergy_corr);
		histos_["sumenergy_wheel0"]->Fill(sumsimhits*1000*cos(angle_rad));
		
		if(hit_maxenergy_corr >= 1.3 && hit_maxenergy_corr <= 2.8){
			histos_["geneta_aftercorr"]->Fill(geneta);
		}
		
	} else {
		histos_["hit_maxenergy_outer"]->Fill(hit_maxenergy_corr);
		histos_["sumenergy_outer"]->Fill(sumsimhits*1000*cos(angle_rad));
	}
	
	if(hit_maxenergy_corr >= 4.6 && hit_maxenergy_corr <= 6.0) {
		histos_["part_geneta_corr"]->Fill(geneta);
		histos_["ieta_iphi_hit"]->Fill(maxEta,maxPhi);
	}
	
	if(hit_maxenergy_corr >= 1.4){
		++i;
		etaphiScatter_ge->SetPoint(i,geneta,genphi);
	}
	
	if(hit_maxenergy_corr < 1.4 && hit_maxenergy_corr > 0.){
		++j;
		etaphiScatter_be->SetPoint(i,geneta,genphi);
	}
	
	if(hit_maxenergy_corr < 0.){
		++k;
		etaphiScatter_ne->SetPoint(i,geneta,genphi);
	}


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
Ho_MuGunAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Ho_MuGunAnalyzer::endJob() 
{
	
	HODeadRegions->Write("HODeadRegions");
	HOSiPMRegions->Write("HOSiPMRegions");
	etaphiScatter_ge->Write("etaphiScatter_ge");
	etaphiScatter_be->Write("etaphiScatter_be");
	etaphiScatter_ne->Write("etaphiScatter_ne");
	theOutputFile->Write();
	theOutputFile->Close();
	std::cout << "i:	" << i << std::endl;
	std::cout << "j:	" << j << std::endl;
	std::cout << "k:	" << k << std::endl;
	
}

// ------------ method called when starting to processes a run  ------------
void 
Ho_MuGunAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Ho_MuGunAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Ho_MuGunAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Ho_MuGunAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Ho_MuGunAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ho_MuGunAnalyzer);
