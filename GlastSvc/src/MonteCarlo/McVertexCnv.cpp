// McVertexCnv.cpp: implementation of the McVertexCnv class.
//
//////////////////////////////////////////////////////////////////////
// $Header$

#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ObjectVector.h"

#include "McVertexCnv.h"
#include "instrument/MCTruth.h"

//! Instantiation of a static factory class used by clients to create instances
static CnvFactory<McVertexCnv> s_factory;
const ICnvFactory& McVertexCnvFactory = s_factory;


//! Create the transient representation of an object.

StatusCode McVertexCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject)   
{

    refpObject = m_vlist= new McVertexCol;
    DataObject* mcroot;
    StatusCode sc;
    sc=dataProvider()->findObject("/Event/MC", mcroot);
    sc=dataProvider()->registerObject(mcroot, "McParticleCol", m_plist=new McParticleCol);

    MCTruth* t = MCTruth::instance();
    if( t ) {
        // do this if MCTruth exists. Otherwise assume MCTree will take care of it
       FluxGenerator * f = dynamic_cast<FluxGenerator*>(t->particle()); // first particle
       addParticle(0, f->m_particles);
    }
  
    return StatusCode::SUCCESS;
}

//! Update object from scratch (not implemented here--does it matter?)
StatusCode McVertexCnv::updateObj(IOpaqueAddress* pAddress, DataObject* pObject)   {
    return StatusCode::SUCCESS;
    
}

//! class ID -- note that it uses the McVertex class ID
const CLID& McVertexCnv::classID()    {
    return McVertex::classID(); 
}

//! Standard Constructor
McVertexCnv::McVertexCnv(ISvcLocator* svc) 
: BaseCnv(classID(), svc)
{

  declareObject("/Event/MC/McVertexCol", objType(), "");
}

//! Standard Destructor
McVertexCnv::~McVertexCnv(){ }

void McVertexCnv::addParticle(McVertex* mother, FluxGenerator::ParticleData* p)
{
    // create a new node in the tree
    McParticle* p1 = new McParticle;
    McVertex*   v1 = new McVertex;
    m_vlist->add(v1);
    m_plist->add(p1);

    // a McVertex is really a track segment
    v1->setInitialPosition(p->m_pos);
    if( p->m_kids.size()>0)  v1->setFinalPosition(p->m_kids[0]->m_pos);
    Hep3Vector momentum=p->m_p;
    int id = p->m_idcode;
    double mass=0;
    switch (abs(id)) { //TODO: use the particle propery database
    case 11:   mass = 0.000511; break;
    case 211:  mass = 0.1395675; break;
    case 13:   mass = 0.105658; break;
    case 2212: mass = 0.938272; break;
    };
    double e2 = mass*mass+momentum.mag2();
    v1->setMcParticle(p1); // associated particle will be set up below
    v1->setInitialFourMomentum( HepLorentzVector(momentum,sqrt(e2)) );
    p1->setParticleID(id); 
    p1->setParticleProperty(id);
    
    //TODO: choose among primaryOrigin, daughterOrigin, decayProduct, showerContents, showerBacksplash
    v1->setVertexType(McVertex::daughterOrigin); 

    v1->setMotherMcParticle( mother!=0? mother->mcParticle() : 0);

    p1->setPrimaryParticleFlag(mother==0);
    p1->setMcVertex(v1);
    if( mother !=0) {
        mother->addDaughterMcParticle(p1);
    } else  m_root = v1;  // make root available for display, etc.

    HepLorentzVector final;
    for( FluxGenerator::ParticleData::ChildList::iterator it= p->m_kids.begin(); 
    it !=p->m_kids.end(); ++it){
        addParticle(v1, *it);
//        final +=  (*it)->m_p; 
    }
    v1->setFinalFourMomentum(final);
    return;
}