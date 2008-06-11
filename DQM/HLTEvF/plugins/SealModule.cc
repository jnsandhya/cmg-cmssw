// $Id: SealModule.cc,v 1.7 2008/05/30 09:25:48 bachtis Exp $

#include "DQM/HLTEvF/interface/HltAnalyzer.h"
#include "DQM/HLTEvF/interface/PathTimerService.h"
#include "DQM/HLTEvF/interface/PathTimerInserter.h"
#include "DQM/HLTEvF/interface/HLTMuonDQMSource.h"
#include "DQM/HLTEvF/interface/HLTTauDQMSource.h"
#include "DQM/HLTEvF/interface/HLTMonElectron.h"
#include "DQM/HLTEvF/interface/FourVectorHLT.h"
#include "DQM/HLTEvF/interface/HLTEventInfoClient.h"

using edm::service::PathTimerService;

#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
DEFINE_FWK_SERVICE(PathTimerService);
//DEFINE_FWK_SERVICE_MAKER(PathTimerService,PathTimerServiceMaker);

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PathTimerInserter);
DEFINE_FWK_MODULE(HltAnalyzer);
DEFINE_FWK_MODULE(HLTMuonDQMSource);
DEFINE_FWK_MODULE(HLTTauDQMSource);
DEFINE_FWK_MODULE(HLTMonElectron);
DEFINE_FWK_MODULE(FourVectorHLT);
DEFINE_FWK_MODULE(HLTEventInfoClient);

