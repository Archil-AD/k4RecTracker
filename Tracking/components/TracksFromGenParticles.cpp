#include "Gaudi/Property.h"

// edm4hep
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackMCParticleLinkCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

// marlin
#include <marlinutil/HelixClass_double.h>

// k4FWCore
#include "k4FWCore/Transformer.h"
#include "k4FWCore/DataHandle.h"
// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

#include <string>

/** @class TracksFromGenParticles
 *
 *  Gaudi transformer that builds an edm4hep::TrackCollection out of an edm4hep::MCParticleCollection.
 *  It just builds an helix out of the genParticle position, momentum, charge and user defined z component of the (constant) magnetic field.
 *  From this helix, different edm4hep::TrackStates (AtIP, AtFirstHit, AtLastHit and AtCalorimeter) are defined. #FIXME for now these trackstates are dummy (copies of the same helix parameters)
 *  This is meant to enable technical development needing edm4hep::Track and performance studies where having generator based trackis is a reasonable approximation.
 *  Possible inprovement:
 *    - Retrieve magnetic field from geometry: const DD4hep::Field::MagneticField& magneticField = detector.field(); DD4hep::DDRec::Vector3D field = magneticField.magneticField(point);
 *    - Properly define different trackStates
 *
 *  @author Brieuc Francois
 */

class TracksFromGenParticles : public Gaudi::Algorithm {

  public:
    TracksFromGenParticles(const std::string& name, ISvcLocator* svcLoc);
    StatusCode initialize();
    StatusCode execute(const EventContext&) const;
    StatusCode finalize();

  private:
    /// Handle for calo hits (input collection)
    mutable DataHandle<edm4hep::MCParticleCollection> m_inputMCParticles{"MCParticles", Gaudi::DataHandle::Reader, this};
    mutable DataHandle<edm4hep::SimTrackerHitCollection> m_inputSimTrackerHits{"SimTrackerHit", Gaudi::DataHandle::Reader, this};
    Gaudi::Property<float> m_Bz{this, "Bz", 2., "Z component of the (assumed constant) magnetic field in Tesla."};
    Gaudi::Property<float> m_RadiusAtCalo{this, "RadiusAtCalo", 2172.8, "Inner radius (in mm) of the calorimeter where the TrackState::AtCalorimeter should be defined."};
    mutable DataHandle<edm4hep::TrackMCParticleLinkCollection> m_links{"TracksFromGenParticlesAssociation", Gaudi::DataHandle::Writer, this};
    mutable DataHandle<edm4hep::TrackCollection> m_tracks{"TracksFromGenParticles", Gaudi::DataHandle::Writer, this};
};

TracksFromGenParticles::TracksFromGenParticles(const std::string& name, ISvcLocator* svcLoc) :
Gaudi::Algorithm(name, svcLoc) {
  declareProperty("InputGenParticles", m_inputMCParticles, "input MCParticles");
  declareProperty("InputSimTrackerHits", m_inputSimTrackerHits, "input SimTrackerHits");
  declareProperty("OutputTracks", m_tracks, "Output tracks");
  declareProperty("OutputMCRecoTrackParticleAssociation", m_links, "MCRecoTrackParticleAssociation");
}

StatusCode TracksFromGenParticles::initialize() {
  StatusCode sc = Gaudi::Algorithm::initialize();
  if (sc.isFailure()) return sc;
  return StatusCode::SUCCESS;
}

StatusCode TracksFromGenParticles::execute(const EventContext&) const {
    // Get the input collection with Geant4 hits
    const edm4hep::MCParticleCollection* genParticleColl = m_inputMCParticles.get();
    const edm4hep::SimTrackerHitCollection* simTrackerHitColl = m_inputSimTrackerHits.get();
    auto outputTrackCollection = new edm4hep::TrackCollection();
    auto MCRecoTrackParticleAssociationCollection = new edm4hep::TrackMCParticleLinkCollection();

    int iparticle = 0;
    for (const auto& genParticle : *genParticleColl) {
      debug() << "Particle decayed in tracker: " << genParticle.isDecayedInTracker() << endmsg;
      debug() << genParticle << endmsg;

      debug() <<"  particle "<<iparticle++<<"  PDG: "<< genParticle.getPDG()  << " energy: "<<genParticle.getEnergy()
              << " charge: "<< genParticle.getCharge() << endmsg;

      // consider only charged particles
      if(genParticle.getCharge() == 0) continue;

      // save SimTrackerHits position and momentum
      std::vector<std::array<double,6> > trackHits;
      for (const auto& hit : *simTrackerHitColl) {
        const edm4hep::MCParticle particle = hit.getParticle();
        std::array<double,6> ahit{hit.x(), hit.y(), hit.z(), hit.getMomentum()[0], hit.getMomentum()[1], hit.getMomentum()[2]};
        if(particle.getVertex() == genParticle.getVertex() && particle.getMomentum() == genParticle.getMomentum()) trackHits.push_back(ahit);
      }

      // consider only particles with at least one SimTrackerHit
      if(trackHits.empty()) continue;

      // sort the hits according to rho
      std::sort(trackHits.begin(), trackHits.end(), [](const std::array<double,6> a, const std::array<double,6> b) {
        double rhoA = std::sqrt(a[0]*a[0] + a[1]*a[1]);
        double rhoB = std::sqrt(b[0]*b[0] + b[1]*b[1]);
        return rhoA < rhoB;
      });

      // Building an helix out of MCParticle properties and B field
      auto helixFromGenParticle = HelixClass_double();
      double genParticleVertex[] = {genParticle.getVertex().x, genParticle.getVertex().y, genParticle.getVertex().z};
      double genParticleMomentum[] = {genParticle.getMomentum().x, genParticle.getMomentum().y, genParticle.getMomentum().z};
      helixFromGenParticle.Initialize_VP(genParticleVertex, genParticleMomentum, genParticle.getCharge(), m_Bz);

      // create a track object
      auto trackFromGen = edm4hep::MutableTrack();

      // TrackState at IP
      auto trackState_IP = edm4hep::TrackState {};
      trackState_IP.location = edm4hep::TrackState::AtIP;
      trackState_IP.D0 = helixFromGenParticle.getD0();
      trackState_IP.phi = helixFromGenParticle.getPhi0();
      trackState_IP.omega = helixFromGenParticle.getOmega();
      trackState_IP.Z0 = helixFromGenParticle.getZ0();
      trackState_IP.tanLambda = helixFromGenParticle.getTanLambda();
      trackState_IP.referencePoint = edm4hep::Vector3f((float)genParticleVertex[0],(float)genParticleVertex[1],(float)genParticleVertex[2]);
      trackFromGen.addToTrackStates(trackState_IP);

      // TrackState at First Hit
      auto trackState_AtFirstHit = edm4hep::TrackState {};
      double posAtFirstHit[] = {trackHits.front()[0],trackHits.front()[1], trackHits.front()[2]};
      double momAtFirstHit[] = {trackHits.front()[3],trackHits.front()[4], trackHits.front()[5]};
      // get extrapolated momentum from the helix with ref point at IP
      helixFromGenParticle.getExtrapolatedMomentum(posAtFirstHit,momAtFirstHit);
      // produce new helix at first hit position
      auto helixAtFirstHit = HelixClass_double();
      helixAtFirstHit.Initialize_VP(posAtFirstHit, momAtFirstHit, genParticle.getCharge(), m_Bz);
      // fill the TrackState parameters
      trackState_AtFirstHit.location = edm4hep::TrackState::AtFirstHit;
      trackState_AtFirstHit.D0 = helixAtFirstHit.getD0();
      trackState_AtFirstHit.phi = helixAtFirstHit.getPhi0();
      trackState_AtFirstHit.omega = helixAtFirstHit.getOmega();
      trackState_AtFirstHit.Z0 = helixAtFirstHit.getZ0();
      trackState_AtFirstHit.tanLambda = helixAtFirstHit.getTanLambda();
      trackState_AtFirstHit.referencePoint = edm4hep::Vector3f((float)posAtFirstHit[0],(float)posAtFirstHit[1],(float)posAtFirstHit[2]);
      trackFromGen.addToTrackStates(trackState_AtFirstHit);

      // TrackState at Last Hit
      auto trackState_AtLastHit = edm4hep::TrackState{};
      double posAtLastHit[] = {trackHits.back()[0],trackHits.back()[1], trackHits.back()[2]};
      double momAtLastHit[] = {trackHits.back()[3],trackHits.back()[4], trackHits.back()[5]};
      // get extrapolated momentum from the helix with ref point at first hit
      helixAtFirstHit.getExtrapolatedMomentum(posAtLastHit,momAtLastHit);
      // produce new helix at last hit position
      auto helixAtLastHit = HelixClass_double();
      helixAtLastHit.Initialize_VP(posAtLastHit, momAtLastHit, genParticle.getCharge(), m_Bz);
      // fill the TrackState parameters
      trackState_AtLastHit.location = edm4hep::TrackState::AtLastHit;
      trackState_AtLastHit.D0 = helixAtLastHit.getD0();
      trackState_AtLastHit.phi = helixAtLastHit.getPhi0();
      trackState_AtLastHit.omega = helixAtLastHit.getOmega();
      trackState_AtLastHit.Z0 = helixAtLastHit.getZ0();
      trackState_AtLastHit.tanLambda = helixAtLastHit.getTanLambda();
      trackState_AtLastHit.referencePoint = edm4hep::Vector3f((float)posAtLastHit[0],(float)posAtLastHit[1],(float)posAtLastHit[2]);
      // attach the TrackStat to the track
      trackFromGen.addToTrackStates(trackState_AtLastHit);

      // TrackState at Calorimeter
      auto trackState_AtCalorimeter = edm4hep::TrackState{};
      double pointAtCalorimeter[] = {0.,0.,0.,0.,0.,0.};
      auto time = helixAtLastHit.getPointOnCircle(m_RadiusAtCalo,posAtLastHit,pointAtCalorimeter);
      double posAtCalorimeter[] = {pointAtCalorimeter[0],pointAtCalorimeter[1],pointAtCalorimeter[2]};
      double momAtCalorimeter[] = {0.,0.,0.};
      // get extrapolated momentum from the helix with ref point at last hit
      helixAtLastHit.getExtrapolatedMomentum(posAtCalorimeter,momAtCalorimeter);
      // produce new helix at calorimeter position
      auto helixAtCalorimeter = HelixClass_double();
      helixAtCalorimeter.Initialize_VP(posAtCalorimeter, momAtCalorimeter, genParticle.getCharge(), m_Bz);
      // fill the TrackState parameters
      trackState_AtCalorimeter.location = edm4hep::TrackState::AtCalorimeter;
      trackState_AtCalorimeter.D0 = helixAtCalorimeter.getD0();
      trackState_AtCalorimeter.phi = helixAtCalorimeter.getPhi0();
      trackState_AtCalorimeter.omega = helixAtCalorimeter.getOmega();
      trackState_AtCalorimeter.Z0 = helixAtCalorimeter.getZ0();
      trackState_AtCalorimeter.tanLambda = helixAtCalorimeter.getTanLambda();
      trackState_AtCalorimeter.referencePoint = edm4hep::Vector3f((float)posAtCalorimeter[0],(float)posAtCalorimeter[1],(float)posAtCalorimeter[2]);
      // attach the TrackStat to the track
      trackFromGen.addToTrackStates(trackState_AtCalorimeter);

      //debug() << trackFromGen << endmsg;
      outputTrackCollection->push_back(trackFromGen);

      // Building the association between tracks and genParticles
      auto MCRecoTrackParticleAssociation = edm4hep::MutableTrackMCParticleLink();
      MCRecoTrackParticleAssociation.setFrom(trackFromGen);
      MCRecoTrackParticleAssociation.setTo(genParticle);
      MCRecoTrackParticleAssociationCollection->push_back(MCRecoTrackParticleAssociation);
    }


    // push the outputTrackCollection to event store
    m_tracks.put(outputTrackCollection);
    m_links.put(MCRecoTrackParticleAssociationCollection);

    debug() << "Output tracks collection size: " << outputTrackCollection->size() << endmsg;

    return StatusCode::SUCCESS;
}

StatusCode TracksFromGenParticles::finalize() { return Gaudi::Algorithm::finalize(); }


DECLARE_COMPONENT(TracksFromGenParticles)
