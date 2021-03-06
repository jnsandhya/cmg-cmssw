import PhysicsTools.HeppyCore.framework.config as cfg
import os

#####COMPONENT CREATOR
from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
kreator = ComponentCreator()

### common MC samples
from CMGTools.RootTools.samples.samples_13TeV_RunIISpring15MiniAODv2 import *

### additional MC samples

### 25 ns

### 50 ns

###

## === 741 RelVals ===
ADD_MJ = kreator.makeMCComponent("ADD_MJ", "/RelValADDMonoJet_d3MD3_13/CMSSW_7_4_1-MCRUN2_74_V9_gensim_740pre7-v1/MINIAODSIM", "CMS", ".*root")
TTLep = kreator.makeMCComponent("TTLep", "/RelValTTbarLepton_13/CMSSW_7_4_1-MCRUN2_74_V9_gensim_740pre7-v1/MINIAODSIM", "CMS", ".*root")
TTbar = kreator.makeMCComponent("TTbar", "/RelValTTbar_13/CMSSW_7_4_1-MCRUN2_74_V9_gensim71X-v1/MINIAODSIM", "CMS", ".*root")
RSGravGaGa = kreator.makeMCComponent("RSGravGaGa", "/RelValRSGravitonToGaGa_13TeV/CMSSW_7_4_1-MCRUN2_74_V9_gensim71X-v1/MINIAODSIM", "CMS", ".*root")

RelVals741 = [ ADD_MJ, TTLep, TTbar, RSGravGaGa ] 



#-----------DATA--------------- 
from CMGTools.RootTools.samples.samples_13TeV_DATA2015 import *

json = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt"

### ---- private samples ----
privEGamma2015Afiles = [ f.strip() for f in open("%s/src/CMGTools/MonoXAnalysis/python/samples/privEGamma_2015A_MINIAOD.txt"  % os.environ['CMSSW_BASE'], "r") ]
privDoubleEG2015Afiles = [ f.strip() for f in open("%s/src/CMGTools/MonoXAnalysis/python/samples/privDoubleEG_2015A_MINIAOD.txt"  % os.environ['CMSSW_BASE'], "r") ]
privHLTPhysics2015Afiles = [ f.strip() for f in open("%s/src/CMGTools/MonoXAnalysis/python/samples/privHLTPhysics_2015A_MINIAOD.txt"  % os.environ['CMSSW_BASE'], "r") ]
privReRecoDoubleEG2015Cfiles = [ f.strip() for f in open("%s/src/CMGTools/MonoXAnalysis/python/samples/privReReco_DoubleEG2015C_MINIAOD.txt" % os.environ['CMSSW_BASE'], "r") ]
privExpress2015Dfiles = [ f.strip() for f in open("%s/src/CMGTools/MonoXAnalysis/python/samples/privExpress_Run2015D_MINIAOD.txt" % os.environ['CMSSW_BASE'], "r") ]
def _grep(x,l): return [ i for i in l if x in i ]
privEGamma2015A_0T = kreator.makePrivateDataComponent('EGamma2015A', '/store/group/dpg_ecal/comm_ecal/data13TeV/EGamma/MINIAOD', _grep('2015A', privEGamma2015Afiles), json )
privDoubleEG2015A_0T = kreator.makePrivateDataComponent('DoubleEG2015A', '/store/group/dpg_ecal/comm_ecal/data13TeV/DoubleEG/MINIAOD', _grep('PAT', privDoubleEG2015Afiles), json )
privHLTPhysics2015A_0T = kreator.makePrivateDataComponent('HLTPhysics2015A', '/store/group/dpg_ecal/comm_ecal/data13TeV/HLTPhysics/MINIAOD', _grep('HLTPhysics', privHLTPhysics2015Afiles), json )
privReRecoDoubleEG2015C = kreator.makePrivateDataComponent('DoubleEG2015C', '/store/group/dpg_ecal/comm_ecal/localreco/data2015_zskim_multifits/miniaod//run2015C_dataPS', _grep('2015C', privReRecoDoubleEG2015Cfiles), json )
privExpress2015D = kreator.makePrivateDataComponent('Express2015D', '/store/group/dpg_ecal/comm_ecal/localreco/data2015D_zskim/miniaod/ExpressPhysics/256/843', _grep('Express', privExpress2015Dfiles), json )
privDataSamples_0T = [ privEGamma2015A_0T, privDoubleEG2015A_0T, privHLTPhysics2015A_0T ]

### ----------------------------- summary ----------------------------------------     

### --- data ---
dataSamples = dataSamples_Run2015D_v4 + dataSamples_Run2015B_05Oct + dataSamples_Run2015D_05Oct

### --- mc ---
#mcSamples_monojet_Asymptotic25ns = TTs + SingleTop + [WJetsToLNu, DYJetsToLL_M50, DYJetsToLL_M10to50 ] + ZJetsToNuNuHT + DYJetsM50HT + WJetsToLNuHT + GJetsHT + QCDHT + DiBosons
mcSamples_monojet_Asymptotic25ns = [ TTJets_LO ] + [WJetsToLNu, DYJetsToLL_M50, DYJetsToLL_M10to50] + ZJetsToNuNuHT + DYJetsM50HT + WJetsToLNuHT + DiBosons
mcSamples_monojet = mcSamples_monojet_Asymptotic25ns 
# 50ns not in miniAOD v2
#mcSamples_monojet_Asymptotic50ns = [ TTJets_50ns, TTJets_LO_50ns, DYJetsToLL_M50_50ns, WJetsToLNu_50ns ] + QCDPt_50ns
#mcSamples_monojet = RelVals740 + RelVals741 + mcSamples_monojet_Asymptotic25ns + mcSamples_monojet_Asymptotic50ns


from CMGTools.TTHAnalysis.setup.Efficiencies import *
dataDir = "$CMSSW_BASE/src/CMGTools/TTHAnalysis/data"

#Define splitting
for comp in mcSamples:
    comp.isMC = True
    comp.isData = False
    comp.splitFactor = 250 #  if comp.name in [ "WJets", "DY3JetsM50", "DY4JetsM50","W1Jets","W2Jets","W3Jets","W4Jets","TTJetsHad" ] else 100
    comp.puFileMC=dataDir+"/puProfile_Summer12_53X.root"
    comp.puFileData=dataDir+"/puProfile_Data12.root"
    comp.efficiency = eff2012
for comp in dataSamples:
    comp.splitFactor = 1000
    comp.isMC = False
    comp.isData = True


if __name__ == "__main__":
   import sys
   if "test" in sys.argv:
       from CMGTools.RootTools.samples.ComponentCreator import testSamples
       testSamples(mcSamples)
