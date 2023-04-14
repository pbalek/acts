#!/usr/bin/env python3
import pathlib, acts, acts.examples, acts.examples.itk
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addPythia8,
    addFatras,
    ParticleSelectorConfig,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    SeedingAlgorithm,
    TruthSeedRanges,
    addCKFTracks,
    CKFPerformanceConfig,
    TrackSelectorRanges,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
    addSeedVertexFinding,
)

ttbar_pu200 = False
u = acts.UnitConstants
geo_dir = pathlib.Path("acts-itk")
outputDir = pathlib.Path.cwd() / "itk_output"
# acts.examples.dump_args_calls(locals())  # show acts.examples python binding calls

detector, trackingGeometry, decorators = acts.examples.itk.buildITkGeometry(geo_dir)
field = acts.examples.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(events=10, numThreads=100, outputDir=str(outputDir))

import sys
print('argument list: %s' % str(sys.argv))

numParticles=int(sys.argv[1])
# particle=int(sys.argv[2])
minpT=float(sys.argv[3])
vtxX,vtxY,vtxZ=float(sys.argv[4]),float(sys.argv[5]),float(sys.argv[6])

if sys.argv[2]=="13":
    particle=acts.PdgParticle.eMuon
if sys.argv[2]=="211":
    particle=acts.PdgParticle.ePionPlus
if sys.argv[2]=="2112":
    particle=acts.PdgParticle.eProton
if sys.argv[2]=="-13":
    particle=acts.PdgParticle.eAntiMuon
if sys.argv[2]=="-211":
    particle=acts.PdgParticle.ePionMinus
if sys.argv[2]=="-2112":
    particle=acts.PdgParticle.eAntiProton

if not ttbar_pu200:
    addParticleGun(
        s,
        MomentumConfig(minpT*u.GeV, 10*minpT*u.GeV, transverse=True),
        EtaConfig(-4.0, 4.0, uniform=True),
        # ParticleConfig(10, acts.PdgParticle.ePionPlus, randomizeCharge=True),
        ParticleConfig(numParticles, particle, randomizeCharge=False),
        vtxGen=acts.examples.GaussianVertexGenerator(
            stddev=acts.Vector4(0.01 * u.mm, 0.01 * u.mm, 0.01 * u.mm, 0.0 * u.ns),
            mean=acts.Vector4(vtxX*u.mm, vtxY*u.mm, vtxZ*u.mm, 0),
        ),
        rnd=rnd,
    )
else:
    addPythia8(
        s,
        hardProcess=["Top:qqbar2ttbar=on"],
        npileup=200,
        vtxGen=acts.examples.GaussianVertexGenerator(
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
            mean=acts.Vector4(0, 0, 0, 0),
        ),
        rnd=rnd,
        outputDirRoot=outputDir,
    )

addFatras(
    s,
    trackingGeometry,
    field,
    rnd=rnd,
    preSelectParticles=ParticleSelectorConfig(
        eta=(-4.0, 4.0), pt=(150 * u.MeV, None), removeNeutral=True
    )
    if ttbar_pu200
    else ParticleSelectorConfig(),
    outputDirRoot=outputDir,
)

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=geo_dir / "itk-hgtd/itk-smearing-config.json",
    outputDirRoot=outputDir,
    rnd=rnd,
)

addSeeding(
    s,
    trackingGeometry,
    field,
    TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-4.0, 4.0), nHits=(9, None))
    if ttbar_pu200
    else TruthSeedRanges(),
    seedingAlgorithm=SeedingAlgorithm.Default,
    *acts.examples.itk.itkSeedingAlgConfig("PixelSpacePoints"),
    geoSelectionConfigFile=geo_dir / "itk-hgtd/geoSelection-ITk.json",
    outputDirRoot=outputDir,
)

addCKFTracks(
    s,
    trackingGeometry,
    field,
    CKFPerformanceConfig(ptMin=1.0 * u.GeV if ttbar_pu200 else 0.0, nMeasurementsMin=6),
    TrackSelectorRanges(pt=(1.0 * u.GeV, None), absEta=(None, 4.0)),
    outputDirRoot=outputDir,
)

addAmbiguityResolution(
    s,
    AmbiguityResolutionConfig(maximumSharedHits=3),
    CKFPerformanceConfig(ptMin=1.0 * u.GeV if ttbar_pu200 else 0.0, nMeasurementsMin=6),
    outputDirRoot=outputDir,
)

addVertexFitting(
    s,
    field,
    vertexFinder=VertexFinder.Iterative,
    outputDirRoot=outputDir,
)

addSeedVertexFinding(
    s,
    outputDirRoot=outputDir,
    outputVertices="fittedSeedVertices"
)

s.run()
