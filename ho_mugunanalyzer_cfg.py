import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")

#from TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff import *
from TrackingTools.TrackAssociator.default_cfi import TrackAssociatorParameterBlock

#process.GlobalTag.globaltag = 'START53_V23::All'
process.GlobalTag.globaltag = 'START53_V23::All'

#FileService for histograms
process.TFileService = cms.Service("TFileService",
    fileName=cms.string(
        'mugun_output.root'
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/data/erdogan/output/temp/*.root'
    )
)

process.demo = cms.EDAnalyzer('Ho_MuGunAnalyzer',
	TrackAssociatorParameters = TrackAssociatorParameterBlock.TrackAssociatorParameters
)


process.p = cms.Path(process.demo)
