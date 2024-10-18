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
    TrackSelectorConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
    addSingleSeedVertexFinding,
)

ttbar_pu200 = False
u = acts.UnitConstants
geo_dir = pathlib.Path("acts-itk")
outputDir = pathlib.Path.cwd() / "itk_output_Hough"
# acts.examples.dump_args_calls(locals())  # show acts.examples python binding calls

import sys
print('argument list: %s' % str(sys.argv))

detector, trackingGeometry, decorators = acts.examples.itk.buildITkGeometry(geo_dir)
field = acts.examples.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))
rnd = acts.examples.RandomNumbers(seed=int(sys.argv[1]))

s = acts.examples.Sequencer(events=3000, numThreads=15, outputDir=str(outputDir))

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
if sys.argv[2]=="1000822080":
    particle=acts.PdgParticle.eLead
if sys.argv[2]=="-13":
    particle=acts.PdgParticle.eAntiMuon
if sys.argv[2]=="-211":
    particle=acts.PdgParticle.ePionMinus
if sys.argv[2]=="-2112":
    particle=acts.PdgParticle.eAntiProton

if not ttbar_pu200:
    if particle!=acts.PdgParticle.eLead:
        addParticleGun(
            s,
            MomentumConfig(minpT*u.GeV, 10*minpT*u.GeV, transverse=True),
            EtaConfig(-4.0, +4.0, uniform=True),
            # ParticleConfig(10, acts.PdgParticle.ePionPlus, randomizeCharge=True),
            ParticleConfig(numParticles, particle, randomizeCharge=True),
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(0.01 * u.mm, 0.01 * u.mm, 0.01 * u.mm, 0.0 * u.ns),
                mean=acts.Vector4(vtxX*u.mm, vtxY*u.mm, vtxZ*u.mm, 0),
            ),
            rnd=rnd,
        )
    else:
        ttbar_pu200=True
        addPythia8(
            s,
            npileup=0,
            beam=acts.PdgParticle.eLead,
            cmsEnergy=5.0*u.TeV,
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(0.5 * u.mm, 0.5 * u.mm, 50.0 * u.mm, 0.0 * u.ns),
                mean=acts.Vector4(vtxX*u.mm, vtxY*u.mm, vtxZ*u.mm, 0),
            ),
            rnd=rnd,
            outputDirRoot=outputDir,
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
    rnd,
    ParticleSelectorConfig(
        rho=(0.0 * u.mm, 28.0 * u.mm),
        absZ=(0.0 * u.mm, 1.0 * u.m),
        eta=(-4.0, 4.0),
        pt=(150 * u.MeV, None),
        removeNeutral=True,
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
    TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-4.0, 4.0), nHits=(9, None)) if ttbar_pu200 else TruthSeedRanges(),
    seedingAlgorithm=SeedingAlgorithm.Default,
    *acts.examples.itk.itkSeedingAlgConfig(acts.examples.itk.InputSpacePointsType.PixelSpacePoints),
    geoSelectionConfigFile=geo_dir / "itk-hgtd/geoSelection-ITk.json",
    outputDirRoot=outputDir,
    # outputDirCsv=outputDir,
)

# addCKFTracks(
#     s,
#     trackingGeometry,
#     field,
#     TrackSelectorConfig(
#         pt=(1.0 * u.GeV if ttbar_pu200 else 0.0, None),
#         absEta=(None, 4.0),
#         nMeasurementsMin=6,
#     ),
#     outputDirRoot=outputDir,
# )

# addAmbiguityResolution(
#     s,
#     AmbiguityResolutionConfig(maximumSharedHits=3),
#     outputDirRoot=outputDir,
# )

# addVertexFitting(
#     s,
#     field,
#     vertexFinder=VertexFinder.Iterative,
#     outputDirRoot=outputDir,
# )

# addSingleSeedVertexFinding(
#     s,
#     outputDirRoot=outputDir,
#     outputVertices="fittedSeedVertices1",
#     ecc=1.0,
# )

# addSingleSeedVertexFinding(
#     s,
#     outputDirRoot=outputDir,
#     outputVertices="fittedSeedVertices2",
#     ecc=2.0,
# )

# addSingleSeedVertexFinding(
#     s,
#     outputDirRoot=outputDir,
#     outputVertices="fittedSeedVertices3",
#     ecc=3.0,
# )

# addSingleSeedVertexFinding(
#     s,
#     outputDirRoot=outputDir,
#     outputVertices="fittedSeedVertices5",
#     ecc=5.0,
# )

addSingleSeedVertexFinding(
    s,
    outputDirRoot=outputDir,
    outputVertices="fittedSeedVerticesHough",
    ecc=-1.0,
)

# addSingleSeedVertexFinding(
#     s,
#     outputDirRoot=outputDir,
#     outputVertices="fittedSeedVertices7",
#     ecc=7.0,
# )

# addSingleSeedVertexFinding(
#     s,
#     outputDirRoot=outputDir,
#     outputVertices="fittedSeedVertices10",
#     ecc=10.0,
# )

# addSingleSeedVertexFinding(
#     s,
#     outputDirRoot=outputDir,
#     outputVertices="fittedSeedVertices15",
#     ecc=15.0,
# )

# addSingleSeedVertexFinding(
#     s,
#     outputDirRoot=outputDir,
#     outputVertices="fittedSeedVertices25",
#     ecc=25.0,
# )

s.run()
