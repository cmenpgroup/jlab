import os
env = Environment(ENV = os.environ)
env.ParseConfig('root-config --cflags --glibs')
t = env.Program(target='ctProcess_omega', source=['ctProcess_omega.cc', 'DetectedParticles.cc', 'IntermediateParticles.cc', 'ReconstructedParticles.cc', 'ParticleList.cc', 'EG2Target.cc', 'EG2Cuts.cc', 'ElectronID.cc', 'EC_geometry.cc', 'OmegaMixedEvent.cc', 'Vertex_Corrections.cc', 'PhotonID.cc', 'ChargedPionID.cc', 'HistManager.cc', 'KineReader.cc', 'PartReader.cc', 'RunCounter.cc'])
Default(t)
