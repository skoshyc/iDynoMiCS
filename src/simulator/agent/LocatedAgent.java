/**
 * \package agent
 * \brief Package of utilities that create and manage agents in the simulation
 * and their participation in relevant reactions.
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 */
package simulator.agent;

import idyno.SimTimer;

import java.util.LinkedList;
import java.awt.Color;

import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;
import simulator.*;
import simulator.geometry.CollisionEngine;
import simulator.geometry.ContinuousVector;
import simulator.geometry.Domain;
import simulator.geometry.EuclideanVector;
import simulator.geometry.Quaternion;
//import simulator.geometry.EuclideanVector;
import simulator.geometry.boundaryConditions.AllBC;
import simulator.geometry.boundaryConditions.BoundaryCyclic;
import java.util.Random;
/**
 * \brief Extends ActiveAgent by adding functionality to control agent grid
 * location, agent shoving, agent death and division, and agent movement.
 *  
 * During each global timestep, agent divisions and agent growth lead to many
 * cases where neighbouring agents will overlap. A relaxation algorithm is
 * used to find iteratively the new overlap-minimising steady state
 * configuration of agent locations at the end of each timestep. 
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre
 * for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK)
 * @author Rob Clegg (rjc096@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK)
 */
public abstract class LocatedAgent extends ActiveAgent implements Cloneable 
{
	/**
	 * Temporary store of the new location this cell will move to.
	 */
	protected static ContinuousVector  _newLoc = new ContinuousVector();
	
	/**
	 * Radius of this agent.
	 */
	public Double _radius = 5.0;
	
	/**
	 * Cell radius including any capsules.
	 */
	public Double _totalRadius = 5.0;
	
	/**
	 * Height of this agent.
	 */
	public Double _height = 10.0;
	
	/**
	 * Cell Height including any capsules.
	 */
	public Double _totalHeight = 10.0;
	
	/**
	 * Direction of the hypha
	 */
	
	//public Double _directionx = Math.random();
	
	//public Double _directiony = Math.random();
	
	public static double distance = 0; //stores the latest calculated distance
	public static double[] intersectionPoints[] = {{0,0,0},{0,0,0}};
	/**
	 * @uml.property  name="intersectionPointsV"
	 * @uml.associationEnd  multiplicity="(0 -1)"
	 */
	public static ContinuousVector[] intersectionPointsV = new ContinuousVector[2];
	
	
	
	/**
	 * Random number generator to add or subtract the directions. Used in getLocationHeight()
	 */
	//public static Random random= new Random();
	
	/**
	 * Radius at which this agent will divide.
	 */
	protected Double _myDivRadius = 0.0;
	
	/**
	 * Radius at which agent death is triggered.
	 */
	protected Double _myDeathRadius = 0.0;
	
	/**
	 * Volume of this agent.
	 */
	protected Double _volume = 0.0;
	
	/**
	 * Cell volume including any capsules.
	 */
	protected Double _totalVolume = 0.0;
	
	/**
	 * Agent position - continuous coordinates.
	 */
	protected ContinuousVector _location = new ContinuousVector();
	
	/***
	 * Agent position with height-continuous coordinates
	 */
	protected ContinuousVector locationHeight = new ContinuousVector();
	
	/***
	 * Agent mid point position with height-continuous coordinates
	 */
	
	protected ContinuousVector center1 = new ContinuousVector(); 
	/**
	 * ContinuousVector noting the move that will be applied to the agents position.
	 */
	protected ContinuousVector _movement = new ContinuousVector();
	
	/**
	 * Direction in which this cell divides.
	 */
	protected ContinuousVector _divisionDirection = new ContinuousVector();
	
	/**
	 * List of neighbouring agents in this agent's vicinity.
	 */
	protected LinkedList<LocatedAgent> _myNeighbors = new LinkedList<LocatedAgent>();

	/**
	 * Index of the agent position on the vectorized grid.
	 */
	protected int _agentGridIndex;
	
	/**
	 * Boolean noting whether this agent is interacting with a surface (true)
	 * or not (false).
	 * 
	 */
	protected Boolean _isAttached = false;

	/**
	 * Detachment priority
	 */
	public Double detPriority = 0.0;

	/**
	 * Stores the simulation time since the last division check
	 */
	public Double _timeSinceLastDivisionCheck = Double.MAX_VALUE;

	/**
	 * Distance based probability from a given neighbour (used in HGT).
	 */
	public Double _distProb = 0.0; 								
	
	/**
	 * Cumulative probability as to whether the plasmid will be transferred.
	 */
	public Double _distCumProb = 0.0; 	


	/**
	 * \brief Constructor used to generate progenitor and initialise an object
	 * to store relevant parameters. 
	 */
	public LocatedAgent()
	{
		super();
		_speciesParam = new LocatedParam();
	}
	
	public void randomizeOrientation()
	{
		if (_agentGrid.is3D)
			{
			//_location = new ContinuousVector(0.0,ExtraMath.random.nextDouble() , ExtraMath.random.nextDouble());//random orientation;
			//locationHeight = new ContinuousVector(0.0, ExtraMath.random.nextDouble() , ExtraMath.random.nextDouble()); 
			locationHeight = new ContinuousVector(_location.x+_height ,_location.y+_height , _location.z+_height); //random orientation;
			center1.sendSum(_location, locationHeight);
			center1.times(0.5);}
		else {
			//_location = new ContinuousVector(ExtraMath.random.nextDouble() , 0.0, 0.0);//random orientation;
			//_location = new ContinuousVector(0.5 , 5.0, 0.0);//random orientation;
			//locationHeight = new ContinuousVector(ExtraMath.random.nextDouble() , ExtraMath.random.nextDouble(), 0.0); //random orientation;
			locationHeight = new ContinuousVector(_location.x+_height,_location.y+_height , _location.z); 
			center1.sendSum(_location, locationHeight);
			center1.times(0.5);
		}
		/*System.out.println(_location);
		System.out.println(locationHeight);
		System.out.println(center1);*/
	}
	
	
	
	
	
	
	/**
	 * \brief Creates a daughter Located Agent by cloning this agent and
	 * parameter objects.
	 * 
	 * @throws CloneNotSupportedException Thrown if the agent cannot be cloned.
	 */
	@Override
	@SuppressWarnings("unchecked")
	public Object clone() throws CloneNotSupportedException
	{
		LocatedAgent o = (LocatedAgent) super.clone();
		o._location = (ContinuousVector) this._location.clone();
		o.locationHeight=(ContinuousVector) this.locationHeight.clone();
		o._movement = (ContinuousVector) this._movement.clone();
		o._divisionDirection = (ContinuousVector)
											this._divisionDirection.clone();
		o._myNeighbors = (LinkedList<LocatedAgent>) this._myNeighbors.clone();
		o._agentGridIndex = this._agentGridIndex;
		return o;
	}
	
	/**
	 * \brief Create a new agent in a specified position.
	 * 
	 * @param position	Vector stating where this agent should be located.
	 */
	public void createNewAgent(ContinuousVector position) 
	{
		try 
		{
			// Get a clone of the progenitor.
			LocatedAgent baby = (LocatedAgent) sendNewAgent();
			baby.giveName();
			baby.updateSize();
			
			this._myDivRadius = getDivRadius();
			baby._myDivRadius = getDivRadius();
			baby._myDeathRadius = getDeathRadius();
			
			// Just to avoid to be in the carrier.
			// TODO Rob 13Mar2015: Is this correct?
			position.x += this._totalRadius;
			
			baby.setLocation(position);
			baby.registerBirth();
		} 
		catch (CloneNotSupportedException e) 
		{
			LogFile.writeError(e, "LocatedAgent.createNewAgent()");
		}
	}
	
	protected void updateOrientationVector(double angle)
	{
		//_headLocation = new ContinuousVector(this._location.x, this._location.y + this._radius - this._capsular_radius, this._location.z); 
		//_tailLocation = new ContinuousVector(this._location.x, this._location.y - this._radius + this._capsular_radius, this._location.z);
		
		ContinuousVector[] vo = new ContinuousVector[2];
		vo[0] = this.center1;
		vo[1] = new ContinuousVector(this.center1.x+1,this.center1.y,this.center1.z+1);
		
		ContinuousVector[] v = new ContinuousVector[2]; 
		v[0] = this.center1;
		
		angle = Math.toRadians(90) ;//Math.random(); //in radians
		
		//rotate head from center
		v[1] = this.locationHeight;
		this.locationHeight = RotateVector(angle,v,vo);
		//System.out.println(locationHeight);
		//rotate tail from center
		v[1] = this._location;
		this._location = RotateVector(angle,v,vo);
		
	}
	
	
	
	
	
	/**
	 * \brief Creates an agent of the specified species and notes the grid in
	 * which this is assigned.
	 *
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param xmlMarkUp	A species mark-up within the specified protocol file.
	 */
	@Override
	public void initFromProtocolFile(Simulator aSim, XMLParser xmlMarkUp) 
	{	
		super.initFromProtocolFile(aSim, xmlMarkUp);
		_myDivRadius = getDivRadius();
		_myDeathRadius = getDeathRadius();
	}
	
	/**
	 * \brief Create an agent using information in a previous state or
	 * initialization file.
	 *
	 * Reads in data from the end of the singleAgentData array and then passes
	 * the remaining values onto the super class.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param singleAgentData	Data from the result or initialisation file
	 * that is used to recreate this agent.
	 */
	@Override
	public void initFromResultFile(Simulator aSim, String[] singleAgentData) 
	{
		/*
		 * Find the position to start at.
		 */
		int nValsRead = 5;
		int iDataStart = singleAgentData.length - nValsRead;
		/*
		 * This is necessary for the case when agents in a biofilm
		 * simulation are transferred into a chemostat.
		 */
		if ( Simulator.isChemostat )
			_location.reset();
		else
		{
			Double newAgentX, newAgentY, newAgentZ;
			newAgentX = Double.parseDouble(singleAgentData[iDataStart]);
			newAgentY = Double.parseDouble(singleAgentData[iDataStart+1]);
			if ( _agentGrid.is3D )
				newAgentZ = Double.parseDouble(singleAgentData[iDataStart+2]);
			else
				newAgentZ = 0.0;
			_location.set(newAgentX, newAgentY, newAgentZ);
		}
		/*
		 * Agent size.
		 */
		_radius      = Double.parseDouble(singleAgentData[iDataStart+3]);
		_totalRadius = Double.parseDouble(singleAgentData[iDataStart+4]);
		/*
		 * These are randomly generated.
		 */
		_myDivRadius = getDivRadius();
		_myDeathRadius = getDeathRadius();
		/*
		 * Now go up the hierarchy with the rest of the data.
		 */
		String[] remainingSingleAgentData = new String[iDataStart];
		for (int i=0; i<iDataStart; i++)
			remainingSingleAgentData[i] = singleAgentData[i];
		super.initFromResultFile(aSim, remainingSingleAgentData);
	}

	
	/**
	 * \brief Called at each time step of the simulation to compute cell
	 * growth, update size, and monitor cell death and division.
	 * 
	 * Also determines whether the agent has reached the size at which it must
	 * divide.
	 */
	@Override
	protected void internalStep()
	{
		/*
		 * Compute mass growth over all compartments.
		 */
		grow();
		/*
		 * Apply this mass growth of all compounds on global radius and mass.
		 */
		updateSize();
		/*
		 * Divide if you have to.
		 */
		if ( willDivide() )
			divide();
		/*
		 * Die if you have to.
		 */
		if ( willDie() )
			die(true);
	}

	/**
	 * \brief Update the radius of the agent from the current mass (and then
	 * the volume) of the agent (EPS included).
	 */
	@Override
	public void updateSize() 
	{
		/* 
		 * Update the totalMass field (sum of the particles masses).
		 */
		updateMass();
		/*
		 * Check the mass is positive.
		 */
		if ( _totalMass < 0.0 )
			LogFile.writeLog("Warning: negative mass on agent "+sendName());
		/*
		 * Sum of (particles masses / particles density).
		 */
		updateVolume();
		/*
		 * Compute radius according to the volume.
		 */
		updateRadius();
		/*
		 * Check if by chance the agent is close enough to a support to be
		 * attached.
		 */
		if ( ! Simulator.isChemostat )
			updateAttachment();
	}

	/**
	 * \brief Captures cell division by making a clone of this agent using the
	 * makeKid method.
	 */
	public void divide()
	{
		try
		{
			makeKid();
		}
		catch (CloneNotSupportedException e)
		{
			LogFile.writeError(e, "LocatedAgent.divide()");
		}
	}

	/**
	 * \brief Determines whether or not a cell has reached the radius where
	 * cell division can be triggered.
	 * 
	 * @return	Boolean stating whether cell division should be triggered
	 * (true) or not (false).
	 */
	public boolean willDivide() 
	{
		/*
		 * This ensures that the checks for when to divide don't occur too
		 * often; at most they will occur at the rate of AGENTTIMESTEP.
		 */
		_timeSinceLastDivisionCheck += SimTimer.getCurrentTimeStep();
		if ( _timeSinceLastDivisionCheck < _agentGrid.getAgentTimeStep() )
			return false;
		_timeSinceLastDivisionCheck = 0.0;
		/*
		 * At this point we will actually check whether to divide.
		 */
		return getRadius(false) > this._myDivRadius;
	}

	/**
	 * \brief Determines whether or not a cell has reached the radius where
	 * cell death can be triggered.
	 * 
	 * @return	Boolean stating whether cell death should be triggered (true)
	 * or not (false).
	 */
	public boolean willDie()
	{
		return (_totalMass < 0.0) || (getRadius(false) <= _myDeathRadius);
	}
	
	/**
	 * \brief With it determined that cell division will occur, create a new
	 * agent from the existing one.
	 * 
	 * @throws CloneNotSupportedException Thrown if the agent cannot be cloned.
	 */
	@Override
	public void makeKid() throws CloneNotSupportedException
	{
		/*
		 * Create the new instance.
		 */
		LocatedAgent baby = (LocatedAgent) sendNewAgent();
		/*
		 * These are all generated randomly.
		 */
		this._myDivRadius = getDivRadius();
		baby._myDivRadius = getDivRadius();
		baby._myDeathRadius = getDeathRadius();
		/*
		 * Update the lineage.
		 */
		recordGenealogy(baby);
		/*
		 * Share mass of all compounds between two daughter cells and compute
		 * new size.
		 */
		divideCompounds(baby, getBabyMassFrac());
		/*
		 * In a chemostat, the daughter cells remain with the coordinates of
		 * their progenitor. Otherwise, compute movement to apply to both
		 * cells and apply it.
		 */
		if ( ! Simulator.isChemostat )
		{
			setDivisionDirection(getInteractDistance(baby)/2);
			baby._movement.subtract(_divisionDirection);
			_movement.add(_divisionDirection);
		}
		/*
		 * Now register the agent inside the guilds and the agent grid.
		 */
		baby.registerBirth();
		baby._netVolumeRate = 0.0;
	}

	/**
	 * \brief On agent division, divides the mass between the old and new
	 * agent, at a specified fraction.
	 * 
	 * @param baby	The new agent, which is inheriting mass.
	 * @param babyMassFrac	The fraction of this agents mass that should be
	 * transferred to the new agent.
	 */
	public void divideCompounds(LocatedAgent baby, Double babyMassFrac)
	{
		/*
		 * Apply babyMassFrac.
		 */
		for (int i = 0; i<particleMass.length; i++)
		{
			baby.particleMass[i] *= babyMassFrac;
			this.particleMass[i] *= 1-babyMassFrac;
		}
		/*
		 * Update radius, mass, volumes and growth rates.
		 */
		updateSize();
		baby.updateSize();
		updateGrowthRates();
		baby.updateGrowthRates();
	}

	/**
	 * \brief On agent division, transfers biomass and EPS between the old and
	 * new agent, at a specified ratio.
	 * 
	 * @param baby	The new agent, which is inheriting mass.
	 * @param babyMassFrac	The ratio of the biomass/EPS that should be 
	 * transferred to the new agent.
	 */
	public void transferCompounds(LocatedAgent baby, Double babyMassFrac)
	{
		Double massToTransfer;
		for (int i = 0; i<particleMass.length; i++)
		{
			massToTransfer = this.particleMass[i] * babyMassFrac;
			baby.particleMass[i] += massToTransfer;
			this.particleMass[i] -= massToTransfer;
		}
		/*
		 * Update radius, mass and volumes.
		 */
		updateSize();
		baby.updateSize();
	}
	
	/**
	 * \brief Set the movement vector that states where to put a newly-created
	 * particle.
	 * 
	 * @param distance	Distance between the this agent and the new agent.
	 */
	public void setDivisionDirection(Double distance)
	{
		Double phi, theta;
		phi = ExtraMath.getUniRandAngle();
		theta = ExtraMath.getUniRandAngle();
		_divisionDirection.x = distance * Math.sin(phi) * Math.cos(theta);
		_divisionDirection.y = distance * Math.sin(phi) * Math.sin(theta);
		if ( _agentGrid.is3D )
			_divisionDirection.z = distance * Math.cos(phi);
		else
			_divisionDirection.z = 0.0;
	}

	/* ______________________ SHOVING ___________________________________ */

	/**
	 * \brief Models a mechanical interaction between two located agents.
	 * 
	 * Implemented by extending classes (LocatedAgent).
	 * 
	 * @param MUTUAL	Whether movement is shared between two agents or
	 * applied only to this one. Set in the protocol file.
	 * @return	The move to be applied once the shoving or pull calculations
	 * have been performed.
	 */
	@Override
	public Double interact(boolean MUTUAL)
	{
		move();
		/*
		 * Rebuild your neighbourhood.
		 */
		getPotentialShovers(getInteractDistance());
		for ( LocatedAgent neighbour : _myNeighbors )
			addPushMovement(neighbour, MUTUAL);
		_myNeighbors.clear();
		return move();
	}

	/**
	 * \brief Mutual shoving : The movement by shoving of an agent is calculated based on the cell overlap and added to the agents movement vector.
	 * 
	 * Mutual shoving : The movement by shoving of an agent is calculated based on the cell overlap and added to the agents movement vector. 
	 * Both agents are moved of half the overlapping distance in opposite directions.
	 * 
	 * @param aNeighbour	 Reference to the potentially shoving neighbour
	 * @param isMutual	Whether movement is shared between two agents or applied only to this one
	 * @param gain	Double noting change in position
	 * @return Boolean stating whether shoving is detected (true) or not (false)
	 */
	public void addPushMovement(LocatedAgent aNeighbor, boolean isMutual)
	{
		/*
		 * Cannot push oneself!
		 */
		if ( aNeighbor == this )
			return;
		ContinuousVector diff=new ContinuousVector();
		Double delta;
		if(this.getStringClass().equals("Fungus")||aNeighbor.getStringClass().equals("Fungus")) {
			//Going to check the angle first
			/* verify intersection of capsules */
			ContinuousVector bactMe = new ContinuousVector(); 
			bactMe.sendDiff(locationHeight,_location);//top-bottom
			ContinuousVector bactHim = new ContinuousVector();
			bactHim.sendDiff(aNeighbor.locationHeight,aNeighbor._location);
			
			
	   		double angle = bactMe.angle(bactHim);
	   		// if vector are in opposite direction the algorithm fails.
	   		if (angle < 0)
	   			bactHim.sendDiff(aNeighbor._location,aNeighbor.locationHeight);	
			diff = computeDifferenceAxis(aNeighbor);
			//System.out.println(diff);
			/*
			 * Compute effective cell-cell distance.
			 */
			 //delta = diff.norm() - (_totalRadius+aNeighbor._totalRadius);
			 //delta =diff.checkSign(diff);	
			delta = diff.norm() - getInteractDistance(aNeighbor);
			if ( delta < 0.0 )
			{
				diff.normalizeVector(delta);
				if ( isMutual )
				{
					diff.times(0.5);
					aNeighbor._movement.add(diff);
					//aNeighbor._location.add(5.0*aNeighbor._movement.x,5.0*aNeighbor._movement.y,5.0*aNeighbor._movement.z);
					//aNeighbor.locationHeight.add(5.0*aNeighbor._movement.x,5.0*aNeighbor._movement.y,5.0*aNeighbor._movement.z);
				}
				this._movement.subtract(diff);
				//this._location.add(5.0*this._movement.x,5.0*this._movement.y,5.0*this._movement.z);
				//this.locationHeight.add(5.0*this._movement.x,5.0*this._movement.y,5.0*this._movement.z);
			}
			 /*if ( delta < 0.0 )
				{
				//calculate translation and rotation
					EuclideanVector forceMe = new EuclideanVector(intersectionPointsV[0],intersectionPointsV[1]);
					double newMag = forceMe.magnitude - (2* _radius);
					forceMe = forceMe.Normalize();
					forceMe = new EuclideanVector(forceMe.start,forceMe.mag_x * newMag, 
					forceMe.mag_y * newMag, forceMe.mag_z * newMag);
					
					//forceMe.Times(0.5f);
						
					//double[] _center = {_location.x,_location.y,_location.z};
					double[] _center = {center1.x,center1.y,center1.z};
					EuclideanVector N = new  EuclideanVector(forceMe.end,_center);
					EuclideanVector T = forceMe.CrossProduct(N);;
					this.rotationAngle += = CollisionEngine.applyForceToCapsule(
					this.center1, new EuclideanVector(_location,locationHeight),
					_radius, forceMe, -1, null);
					torque = T;torque.Plus(T);
					//System.out.println(this.rotationAngle+"");
					double gain=0.1; //from iDynoBacillus.simulator.agent.LocatedAgent.addSpringAttachment			
					if (isMutual) {
						//forceMe.Times(0.5f);
						this.rotationAngle *= 0.5;
						forceMe = forceMe.Times(gain);
						this._movement.add(forceMe.mag_x,forceMe.mag_y,forceMe.mag_z);
						//System.out.println(this._movement);
						
						aNeighbor.rotationAngle -= rotationAngle; //= rotationAngle;
						aNeighbor.torque = T;aNeighbor.torque.Minus(T);

						aNeighbor.rotationAngle *= 0.5;
						aNeighbor._movement.subtract(forceMe.mag_x,forceMe.mag_y,forceMe.mag_z);
						//System.out.println(aNeighbor._movement);
						
					} else {
						this.rotationAngle *= 0.5;
						forceMe = forceMe.Times(gain);
						this._movement.add(forceMe.mag_x,forceMe.mag_y,forceMe.mag_z);
						//this._location.add(_movement.x,_movement.y,_movement.z);
						//this.locationHeight.add(_movement.x,_movement.y,_movement.z);
					}
				}*/
		}
		else {
		
		/*
		 * Find the vector from your neighbour's cell centre to your cell
		 * centre. This is for spherical agents
		 */
			diff = computeDifferenceVector(aNeighbor);
		/*
		 * Compute effective cell-cell distance.
		 */
			delta = diff.norm() - getInteractDistance(aNeighbor);
			
		
		/*System.out.println(aNeighbor._movement);
		System.out.println(this._movement);*/
		/*
		 * Apply the shoving calculated. If it's mutual, apply half to each.
		 */
		if ( delta < 0.0 )
		{
			diff.normalizeVector(delta);
			if ( isMutual )
			{
				diff.times(0.5);
				aNeighbor._movement.add(diff);
			}
			this._movement.subtract(diff);
		}
	}
		/*System.out.println(diff);
		System.out.println(delta);*/
	}

	/**
	 * \brief Pulling : The movement of agents by a shrinking biofilm. Move calculated and added to the agents movement vector.
	 * 
	 * The movement of agents by a shrinking biofilm. Move calculated and added to the agents movement vector. 
	 * 
	 * TODO Not currently used... consider deleting?
	 * 
	 * @param aNeighbor	 Reference to the potentially shoving neighbour
	 * @param isMutual	Whether movement is shared between two agents or applied only to this one
	 * @param gain	Double noting change in position
	 * @return Boolean stating whether pulling is detected (true) or not (false)
	 */
	public void addSpringMovement(LocatedAgent aNeighbor, boolean isMutual)
	{
		/*
		 * Cannot push oneself!
		 */
		if ( aNeighbor == this )
			return;
		/*
		 * Build the escape vector and find the distance between you and your
		 * neighbour. 
		 */
		ContinuousVector diff = computeDifferenceVector(aNeighbor);
		/*
		 * Compute effective cell-cell distance. This part differs from 
		 * addPushMovement() in that the 
		 */
		Double delta = diff.norm() - getInteractDistance(aNeighbor);
		if (delta > _totalRadius)
			return;
		// TODO Rob 13Mar2015: where does this 5 come from?
		if ( delta > 0.0 )
			delta *= Math.exp(-delta * 5 / _totalRadius);
		/*
		 * Apply the shoving calculated. If it's mutual, apply half to each.
		 */
		diff.normalizeVector(delta);
		if ( isMutual )
		{
			diff.times(0.5);
			aNeighbor._movement.add(diff);
		} 
		this._movement.subtract(diff);
	}

	/**
	 * \brief Computes the shortest vector between this agent and a position
	 * given as a ContinuousVector. Assumes cyclic boundaries.
	 * 
	 * If the vector is all zero's, returns a vector of random direction and
	 * length = 0.01 * radius.
	 * 
	 * TODO Can we do this without assuming cyclic boundaries? I.e. actually
	 * check..
	 * 
	 * @param position	ContinuousVector of position to calculate distance to.
	 * @return The shortest movement vector to go from a to b, taking into
	 * account the cyclic boundary.
	 * @see addOverlapMovement
	 * @see addPullMovement works in 2 and 3D
	 */
	public ContinuousVector computeDifferenceVector(ContinuousVector position)
	{
		Double gridLength;
		ContinuousVector diff = new ContinuousVector();
		diff.sendDiff(_location, position);
		/*
		 * Check periodicity in X.
		 */
		gridLength = _species.domain.length_X;
		if ( Math.abs(diff.x) > 0.5 * gridLength )
			diff.x -= Math.signum(diff.x) * gridLength;
		/*
		 * Check periodicity in Y.
		 */
		gridLength = _species.domain.length_Y;
		if ( Math.abs(diff.y) > 0.5 * gridLength )
			diff.y -= Math.signum(diff.y) * gridLength;
		/*
		 * Check periodicity in Z.
		 */
		if (_agentGrid.is3D)
		{
			gridLength = _species.domain.length_Z;
			if (Math.abs(diff.z) > 0.5 * gridLength)
				diff.z -= Math.signum(diff.z) * gridLength;
		}
		/*
		 * If this is a zero vector, give it random direction and a norm of
		 * 0.01 * radius.
		 */
		if ( diff.isZero() )
		{
			diff.alea(_agentGrid.is3D);
			diff.normalizeVector(0.01*_radius);
		}
		return diff;
	}
	
	
	/**
	 * \brief Computes the shortest vector between two axes of the cylinders. Assumes cyclic boundaries.
	 * The algorithm is based on "On Fast Computation of distance
	 * between line segments" by Lumelsky,1985.
	 * If the vector is all zero's, returns a vector of random direction and
	 * length = 0.01 * radius.
	 * 
	 
	 * @params LocatedAgent aLoc1, aLoc2
	 * @return The shortest movement vector to go from axis a to axis b, taking into
	 * account the cyclic boundary.
	 **/
	public ContinuousVector computeDifferenceAxis(LocatedAgent aLoc)
	{
		Double gridLength;
		ContinuousVector diff = new ContinuousVector();
		ContinuousVector d1 = new ContinuousVector();
		ContinuousVector d2 = new ContinuousVector();
		ContinuousVector d12 = new ContinuousVector();
		ContinuousVector c1 = new ContinuousVector();
		ContinuousVector c2 = new ContinuousVector();
		Double t;
		Double u;
		
		
		d1.sendDiff(locationHeight,_location);
		d2.sendDiff(aLoc.locationHeight,aLoc._location);
		d12.sendDiff(aLoc._location,_location);
		
		Double R = d1.prodScalar(d2);
		Double D1 = d1.prodScalar(d1);
		Double D2 = d2.prodScalar(d2);
		Double S1 = d1.prodScalar(d12);
		Double S2 = d2.prodScalar(d12);
		if(D1<=0 && D2 <=0){
			t=0.0;
			u=0.0;
			diff.sendDiff(_location,aLoc._location);
			double[] C1= {_location.x,_location.y,_location.z};
			double[] C2= {aLoc._location.x,aLoc._location.y,aLoc._location.z};
			/**/
			intersectionPoints[0] = C1;
			intersectionPoints[1] = C2;
			intersectionPointsV[0] = new ContinuousVector(C1[0],C1[1],C1[2]);
			intersectionPointsV[1] = new ContinuousVector(C2[0],C2[1],C2[2]);
			/**/
		}
		
		if(D1<=0){ //First line segment is a point
			t=0.0;
			u=-S2/D2;
			u=Clamp(u,0.0,1.0);
		}
		else{
			if(D2<=0){ // second line segment is a point
				u=0.0;
				t=S1/D1;
				t=Clamp(t,0.0,1.0);
			}
			else{
				//Both line segments are not points
				Double denom = (D1*D2)-(R*R);
				if(denom != 0){ 
					t=(S1*D2-S2*R)/denom;
					t=Clamp(t,0.0,1.0);
				}
				else{ // Lines are parallel
					t=0.0;
				}
				u=(t*R-S2)/D2;
				//To check if u is in [0,1]
				if(u<0.0){
					u=0.0;
					t=S1/D1;
					t=Clamp(t,0.0,1.0);
				}
				else{
					if(u>1.0){
						u=1.0;
						t=(R+S1)/D1;
						t=Clamp(t,0.0,1.0);
					}
				}
				
			}
			d1.times(t);
			c1.sendSum(_location, d1);
			d2.times(u);
			c2.sendSum(aLoc._location, d2);
			diff.sendDiff(c1,c2);
			double[] C1= {locationHeight.x*t,locationHeight.y*t,locationHeight.z*t}; //p1+d1*s;
			double[] C2= {aLoc.locationHeight.x*t,aLoc.locationHeight.y*t,aLoc.locationHeight.z*t}; //p2+d2*t;
			
			/* extras */
			intersectionPoints[0] = C1;
			intersectionPoints[1] = C2;
			intersectionPointsV[0] = new ContinuousVector(C1[0],C1[1],C1[2]);
			intersectionPointsV[1] = new ContinuousVector(C2[0],C2[1],C2[2]);
			return diff;
		}
		/* * Check periodicity in X.*/
		 
		gridLength = _species.domain.length_X;
		if ( Math.abs(diff.x) > 0.5 * gridLength )
			diff.x -= Math.signum(diff.x) * gridLength;
		
		 /** Check periodicity in Y.*/
		 
		gridLength = _species.domain.length_Y;
		if ( Math.abs(diff.y) > 0.5 * gridLength )
			diff.y -= Math.signum(diff.y) * gridLength;
		
		/* * Check periodicity in Z.*/
		 
		if (_agentGrid.is3D)
		{
			gridLength = _species.domain.length_Z;
			if (Math.abs(diff.z) > 0.5 * gridLength)
				diff.z -= Math.signum(diff.z) * gridLength;
		}
		
		/* * If this is a zero vector, give it random direction and a norm of
		 * 0.01 * radius.*/
		 
		if ( diff.isZero() )
		{
			diff.alea(_agentGrid.is3D);
			diff.normalizeVector(0.01*_radius);
		}
		return diff;
		
	}
	
	/* 
	 * Clamp function relates to equation 12 in the Lumelsky paper
	 * Restricts t and u to lie within [min,max]
	 * param: variable, min, max
	 * returns variable restricted to the interval
	 */
	public Double Clamp(Double var,Double min,Double max) {
		if(var < min) {var=min;}
		if(var > max) {var=max;}
		return var;
	}
	
	
	/**
	 * 
	 * @param aLoc
	 * @return
	 */
	public ContinuousVector computeDifferenceVector(LocatedAgent aLoc)
	{
		/*if(aLoc.getStringClass().equals("Fungus")) {
			return computeDifferenceAxis(aLoc._location,aLoc.locationHeight);
		}
		else {}*/
		return computeDifferenceVector(aLoc._location);
	}
	
	
	
	
	
	/**
	 * \brief Find neighbouring agents in a range around you.
	 * 
	 * @param radius	The distance to search around the agent location.
	 */
	public void getPotentialShovers(Double radius)
	{
		_agentGrid.getPotentialShovers(_agentGridIndex, radius, _myNeighbors);
	}

	/**
	 * \brief Pick a random neighbour from the _myNeigbors collection.
	 * 
	 * If used multiple times without changing _myNeighbours, this will be
	 * random selection WITH repetition.
	 * 
	 * @return	A randomly picked neighbour (LocatedAgent object) from the
	 * list of neighbours.
	 */
	public LocatedAgent pickNeighbor()
	{
		if ( _myNeighbors.isEmpty() )
			return null;
		return _myNeighbors.get(ExtraMath.getUniRandInt(_myNeighbors.size()));
	}

	/**
	 * \brief Find siblings of this agent in the immediate surroundings.
	 * 
	 * @param indexSpecies	The index used to reference this species in the
	 * simulation dictionary.
	 */
	public void findCloseSiblings(int indexSpecies) 
	{
		Double shoveDist;
		LocatedAgent aNb;
		/*
		 * Find and count neighbours.
		 */
		getPotentialShovers(getInteractDistance());
		int nNb = _myNeighbors.size();
		/*
		 * Loop through them, only re-appending them to the neighbour list
		 * if they are: (1) different to this agent, (2) the same species as 
		 * this agent, and (3) close enough to this agent.  
		 */
		for ( int iNb = 0; iNb < nNb; iNb++ )
		{
			aNb = _myNeighbors.removeFirst();
			if ( aNb == this || indexSpecies != aNb.speciesIndex)
				continue;
			shoveDist = 2 * (getShoveRadius() + aNb.getShoveRadius());
			if ( getDistance(aNb) <= shoveDist )
				_myNeighbors.addLast(aNb);
		}
	}
	//double rotationAngle = 0;
	double rotationAngle = Math.toRadians(Math.random() * 360);
    
	/**
	 * @uml.property  name="torque"
	 * @uml.associationEnd  
	 */
	EuclideanVector torque = new EuclideanVector(_location,_location);
	
	/**
	 * Apply the rotation angle stored 
	 */
	public void rotate()
	{
		//System.out.println(rotationAngle);
		if (rotationAngle > 0 && torque.magnitude > 0 && 
				rotationAngle < 0.6)
		{ 
			ContinuousVector[] center_head = {center1,locationHeight};
			ContinuousVector[] center_tail = {center1,_location};
			ContinuousVector[] T1 = {
	        		new ContinuousVector(torque.start[0],torque.start[1],torque.start[2]),
	        		new ContinuousVector(torque.start[0]+torque.mag_x,
	        				torque.start[1] + torque.mag_y,
	        				torque.start[2] + torque.mag_z)};
	     
		    ContinuousVector newHead = RotateVector(rotationAngle,center_head,T1);
		    ContinuousVector newTail = RotateVector(rotationAngle,center_tail,T1);
	   
		    
		    if (!Double.isNaN(newHead.x) && !Double.isNaN(newTail.x))
		    {
		    	this.locationHeight = newHead; 
		    	this._location = newTail;
		    }
		}
		
		//rotationAngle = 0;
		//torque = new EuclideanVector(_location,_location);
	}
	
	/**
	 * \brief With the agent move calculated, apply this movement, taking care
	 * to respect boundary conditions.
	 * 
	 * @return Distance moved relative to total radius.
	 */
	@Override
	public Double move()
	{
		/*
		 * Check the movement is valid.
		 */
		if ( ! _movement.isValid() )
		{
			LogFile.writeLog("Incorrect movement coordinates");
			_movement.reset();
		}
		/*
		 * Check we're not trying to move in the Z direction in 2D.
		 */
		if ( !(_agentGrid.is3D) && !(_movement.z.equals(0.0)) )
		{
			_movement.z = 0.0;
			_movement.reset();
			LogFile.writeLog("Agent tried to move in Z direction!");
		}
		/*
		 * No movement planned, finish here.
		 */
		if (_movement.isZero())
			return 0.0;
		
		//Check to see if the boundaries are crossed
		checkBoundariesTailHead();
		_location.add(5.0*_movement.x,5.0*_movement.y,5.0*_movement.z);
		locationHeight.add(5.0*_movement.x,5.0*_movement.y,5.0*_movement.z);
		
		/*_location.add(_movement.x,_movement.y,_movement.z);
		locationHeight.add(_movement.x,_movement.y,_movement.z);*/
		
		/* Now apply the movement.*/
		 
		
		setLocation(getVerifiedLocationFromMovement(_movement));
		if(_species.speciesClass.equals("Fungus")) {
			setLocationHeight(getVerifiedLocationHeightFromMovement(_movement));
			}
		_agentGrid.registerMove(this);
		/*
		 * Calculate how far we've traveled relative to the total radius.
		 */
		Double delta = _movement.norm();
		_movement.reset();
		if(_species.speciesClass.equals("Fungus")) {
			//System.out.println(delta/_height);
			return delta/_height;
		}
		else {return delta/_totalRadius;}
		//return delta/_totalRadius;
	}

	
	public void checkBoundariesTailHead() {
		
		ContinuousVector _newTailLoc = new ContinuousVector();
		_newTailLoc.set(_location);
		_newTailLoc.add(_movement);
		
		ContinuousVector _newHeadLoc = new ContinuousVector();
		_newHeadLoc.set(locationHeight);
		_newHeadLoc.add(_movement);
		
		AllBC aBoundaryT = getDomain().testCrossedBoundary(_newTailLoc);
		AllBC aBoundaryH = getDomain().testCrossedBoundary(_newHeadLoc);
		
		boolean testHead = (aBoundaryH!=null);
		boolean testTail = (aBoundaryT!=null);
		
		
		if (testHead)
		{
			 EuclideanVector force = 
				 new EuclideanVector(locationHeight,aBoundaryH.getOrthoProj(_location));

			 //we need a vector pointing inside with the size of the capsular radius
			 ContinuousVector forceNormal = new ContinuousVector(0.0,0.0,0.0); 
			 //In the iDynoBacillus code the getNormalInside() fn has a cont. vector as input
			 forceNormal = aBoundaryH.getShape().getNormalInside();
			 forceNormal.normalizeVector();
			 forceNormal.times(_radius);
			 force.mag_x += forceNormal.x;
			 force.mag_y += forceNormal.y;
			 force.mag_z += forceNormal.z;
			 
			 _movement.add(force.getContinuousVector());
			 
			
			 double[] _center = {center1.x,center1.y,center1.z};
				EuclideanVector N = new  EuclideanVector(force.end,_center);
				EuclideanVector T = force.CrossProduct(N);;
				this.rotationAngle += CollisionEngine.applyForceToCapsule(
						this.center1, new EuclideanVector(_location,locationHeight),
						_radius, force, -1, null);
				torque = torque.Plus(T);
		}
		
		if (testTail)
		{
			 EuclideanVector force = 
				 new EuclideanVector(_location,aBoundaryT.getOrthoProj(_location));
		
		
			 //we need a vector pointing inside with the size of the capsular radius
			 ContinuousVector forceNormal = new ContinuousVector(0.0,0.0,0.0); 
			 forceNormal = aBoundaryT.getShape().getNormalInside();
			 forceNormal.normalizeVector();
			 forceNormal.times(_radius);
			 force.mag_x += forceNormal.x;
			 force.mag_y += forceNormal.y;
			 force.mag_z += forceNormal.z;
			 //force.Plus(_capsular_radius,_capsular_radius,_capsular_radius);
			 
			 _movement.add(force.getContinuousVector());
			 //_location.add(force.getContinuousVector());
			 
			 double[] _center = {center1.x,center1.y,center1.z};
				EuclideanVector N = new  EuclideanVector(force.end,_center);
				EuclideanVector T = force.CrossProduct(N);;
				this.rotationAngle += CollisionEngine.applyForceToCapsule(
						this.center1, new EuclideanVector(_location,locationHeight),
						_radius, force, -1, null);
				torque = torque.Plus(T);
		}

	}
	
	
	/**
	 * \brief Used by the move method to determine if an agent's move crosses
	 * any of the domain's boundaries.
	 */
	public ContinuousVector getVerifiedLocationFromMovement(ContinuousVector movement) {
		// Search a boundary which will be crossed
		ContinuousVector newLoc = new ContinuousVector(_location);
		newLoc.add(movement);
		return getVerifiedLocation(newLoc);
	}
	
	public ContinuousVector getVerifiedLocationHeightFromMovement(ContinuousVector movement) {
		// Search a boundary which will be crossed
		ContinuousVector newLoc = new ContinuousVector(locationHeight);
		newLoc.add(movement);
		return getVerifiedLocation(newLoc);
	}
	
	public ContinuousVector getVerifiedLocation(ContinuousVector location) {
		AllBC aBoundary = getDomain().testCrossedBoundary(getRadius(true),location);
		int nDim = (_agentGrid.is3D ? 3 : 2);
		int counter = 0;

		/*
		 * Test all boundaries and apply corrections according to crossed
		 * boundaries.
		 */
		while (aBoundary != null || (counter > nDim))
		{
			counter++;
			aBoundary.applyBoundary(this, location);
			aBoundary = getDomain().testCrossedBoundary(getRadius(true),location);
		}
		return location;
	}

	
	/**
	 * \brief Add the reacting concentration of an agent to the received grid
	 * 
	 * Add the reacting concentration of an agent to the received grid
	 * 
	 * @param aSpG	Spatial grid used to sum catalysing mass
	 * @param catalystIndex	Index of the compartment of the cell supporting the reaction
	 */
	@Override
	public void fitMassOnGrid(SpatialGrid aSpG, int catalystIndex)
	{
		if (isDead)
			return;

		Double value = particleMass[catalystIndex]/aSpG.getVoxelVolume();
		if ( ! Double.isFinite(value) )
			value = 0.0;
		aSpG.addValueAt(value, _location);
	}

	/**
	 * \brief Add the total concentration of an agent on received grid
	 * 
	 * Add the total concentration of an agent on received grid
	 * 
	 * @param aSpG	Spatial grid used to sum total mass
	 */
	@Override
	public void fitMassOnGrid(SpatialGrid aSpG) 
	{
		if (isDead)
			return;

		Double value = _totalMass/aSpG.getVoxelVolume();
		if ( ! Double.isFinite(value) )
			value = 0.0;
		aSpG.addValueAt(value, _location);
	}

	/**
	 * \brief Add the total volume rate of an agent on received grid
	 * 
	 * Add the total volume rate of an agent on received grid
	 * 
	 * @param aSpG	Spatial grid used to sum volume
	 */
	public void fitVolRateOnGrid(SpatialGrid aSpG)
	{
		Double value = _netVolumeRate/aSpG.getVoxelVolume();
		if ( ! Double.isFinite(value) )
			value = 0.0;
		try
		{
			aSpG.addValueAt(value, _location);
		}
		catch (ArrayIndexOutOfBoundsException e)
		{
			LogFile.writeLogAlways("Could not put LocatedAgent mass on grid");
			LogFile.writeLogAlways("Problem with location "
													+_location.toString());
			System.exit(-1);
		}
	}

	/**
	 * \brief Add the reaction/growth rate of an agent on received grid, for a specified reaction
	 * 
	 * Add the total reaction/growth rate of an agent on received grid, for a specified reaction
	 * 
	 * @param aRateGrid	Spatial grid used to store total reaction rate
	 * @param reactionIndex	Index of this declared reaction in the simulation dictionary
	 */
	@Override
	public void fitReacRateOnGrid(SpatialGrid aRateGrid, int reactionIndex)
	{
		if (isDead)
			return;
		
		// growthRate is in [fgX.hr-1] so convert to concentration:
		// [fgX.um-3.hr-1 = gX.L-1.hr-1]
		Double value = growthRate[reactionIndex]/aRateGrid.getVoxelVolume();

		if ( ! Double.isFinite(value) )
			value = 0.0;

		aRateGrid.addValueAt(value, _location);
	}

	/* _______________ FILE OUTPUT _____________________ */

	/**
	 * \brief Used in creation of results files - specifies the header of the columns of output information for this agent
	 * 
	 * Used in creation of results files - specifies the header of the columns of output information for this agent
	 * 
	 * @return	String specifying the header of each column of results associated with this agent
	 */
	@Override
	public StringBuffer sendHeader()
	{
		// return the header file for this agent's values after sending those for super
		StringBuffer tempString = super.sendHeader();
		
		// location info and radius
		if(_species.speciesClass.equals("Fungus")) {
		tempString.append(",locationX,locationY,locationZ,locationHeightX,locationHeightY,locationHeightZ,radius,totalRadius,height,totalHeight");}
		else {
			tempString.append(",locationX,locationY,locationZ,radius,totalRadius");
		}
		
		return tempString;
	}

	/**
	 * \brief Used in creation of results files - creates an output string of information generated on this particular agent
	 * 
	 * Used in creation of results files - creates an output string of information generated on this particular agent
	 * 
	 * @return	String containing results associated with this agent
	 */
	@Override
	public StringBuffer writeOutput()
	{
		// write the data matching the header file
		StringBuffer tempString = super.writeOutput();
		//ContinuousVector locht=getLocationHeight();
		//System.out.println(locationHeight);
		if(_species.speciesClass.equals("Fungus")) {
		// location info and radius
		tempString.append(","+_location.x+","+_location.y+","+_location.z+","+locationHeight.x+","+locationHeight.y+","+locationHeight.z+",");
		tempString.append(_radius+","+_totalRadius+","+_height+","+_totalHeight);}
		else {
			// location info and radius
			tempString.append(","+_location.x+","+_location.y+","+_location.z+",");
			tempString.append(_radius+","+_totalRadius);
		}
		
		return tempString;
	}

	/* _______________ RADIUS, MASS AND VOLUME _____________________ */

	/**
	 * \brief Compute the volume on the basis of the mass and density of different compounds defined in the cell
	 * 
	 * Compute the volume on the basis of the mass and density of different compounds defined in the cell
	 */
	public void updateVolume()
	{
		_volume = 0.0;
		for (int i = 0; i<particleMass.length; i++) {
			_volume += particleMass[i]/getSpeciesParam().particleDensity[i];
		}
		_totalVolume = _volume;
	}

	/**
	 * \brief Compute the radius on the basis of the volume The radius evolution is stored in deltaRadius (used for shrinking)
	 * 
	 * Compute the radius on the basis of the volume The radius evolution is stored in deltaRadius (used for shrinking)
	 */
	public void updateRadius() {

		//sonia:chemostat 22.02.2010
		if(Simulator.isChemostat || _species.domain.is3D){
			_radius = ExtraMath.radiusOfASphere(_volume);
			_totalRadius = ExtraMath.radiusOfASphere(_totalVolume);
		}else{
			
			if(_species.speciesClass.equals("Fungus")) {
				_height = ExtraMath.lengthOfACylinder(_volume,
						_radius);
				_totalHeight = ExtraMath.lengthOfACylinder(_totalVolume,
						_totalRadius);
				double magnitude = locationHeight.distance(_location);
				if (magnitude == 0)
				{
					randomizeOrientation();
					 magnitude = locationHeight.distance(_location);
				}
				double magX = ((locationHeight.x - _location.x) / magnitude) * (_height) ;
				double magY = ((locationHeight.y - _location.y) / magnitude) * (_height) ;
				double magZ = ((locationHeight.z - _location.z) / magnitude) * (_height) ;
				
				Double xlocheight = _location.x+magX;
				//Double ylocheight =_location.y;
				//Double xlocheight = _location.x;
				Double ylocheight =_location.y+magY;
				Double zlocheight = _location.z+magZ;
				locationHeight.set(xlocheight, ylocheight, zlocheight);
				
				
			}else {
				_radius = ExtraMath.radiusOfACylinder(_volume,
				_species.domain.length_Z);
				_totalRadius = ExtraMath.radiusOfACylinder(_totalVolume,
				_species.domain.length_Z);
			}
			
		}
		/*System.out.println("\t Volume, radius, height \n");
		System.out.println(_volume);
		System.out.println(_radius);
		System.out.println(_height);
		*/
	}

	/**
	 * \brief Update the attachment, checking if this agent is close enough to
	 * any boundaries.
	 * 
	 * TODO Rob 13Mar2015: Where does this 3 come from?!
	 * 
	 * @return	Boundary that has been crossed.
	 */
	public void updateAttachment()
	{
		for (AllBC aBoundary : getDomain().getAllBoundaries())
			if ( aBoundary.isSupport() &&
						aBoundary.getDistance(_location) <= 3 * _totalRadius )
			{
				_isAttached = true;
				return;
			}
	}
	
	/**
	 * \brief Add movement to the ContinuousVector storing the agents move.
	 * 
	 * @param aMove	ContinuousVector to add to the movement vector.
	 */
	public void addMovement(ContinuousVector aMove)
	{
		this._movement.add(aMove);
	}
	
	/**
	 * \brief Return the set of parameters associated with this agent
	 * (LocatedParam object).
	 * 
	 * @return LocatedParam object of parameters associated with this agent.
	 */
	@Override
	public LocatedParam getSpeciesParam()
	{
		return (LocatedParam) _speciesParam;
	}
	
	/**
	 * \brief Return the volume of this agent, with or without the capsule.
	 *  
	 * @param withCapsule	Boolean noting whether any capsule should be
	 * included in this calculation.
	 * @return	Double specifying the volume of this agent.
	 */
	public Double getVolume(boolean withCapsule)
	{
		return withCapsule ? _totalVolume : _volume;
	}
	
	/**
	 * \brief Return the radius of this agent, with or without the capsule.
	 * 
	 * @param withCapsule	Boolean noting whether any capsule should be
	 * included in this calculation.
	 * @return	Double specifying the radius of this agent
	 */
	public Double getRadius(boolean withCapsule)
	{
		
		if(_species.speciesClass.equals("Fungus")) {
		return (withCapsule ? _totalHeight : _height);
		}
		else {
			return (withCapsule ? _totalRadius : _radius);
		}
	}
	
	/**
	 * \brief Return the mass of this agent, with or without the capsule.
	 *  
	 * @param withCapsule	Boolean noting whether any capsule should be
	 * included in this calculation.
	 * @return	Double specifying the mass of this agent.
	 */
	public Double getMass(Boolean withCapsule)
	{
		return (withCapsule ? _totalMass : _totalMass);
	}
	
	/**
	 * \brief Report whether this cell has any EPS.
	 * 
	 * @return	Boolean noting whether this cell has any EPS.
	 */
	public Boolean hasEPS() 
	{
		return false;
	}
	
	/**
	 * \brief Return the shove factor to be used in shoving for this
	 * species of agent.
	 * 
	 * @return	Double specifying the shove factor that will be applied.
	 */
	public Double getShoveFactor()
	{
		return ((LocatedParam) _speciesParam).shoveFactor;
	}

	/**
	 * \brief Return the shove radius to be used in shoving for this
	 * species of agent.
	 * 
	 * @return	Double specifying the shove radius that will be applied.
	 */
	public Double getShoveRadius()
	{
		if(_species.speciesClass.equals("Fungus")) {
			return _height * getShoveFactor(); 
		}
		else {
			return _totalRadius * getShoveFactor();
			}
		
	}
	
	/**
	 * \brief Return the shove limit to be used in shoving for this
	 * species of agent.
	 * 
	 * @return	Double specifying the shove limit that will be applied.
	 */
	public Double getShoveLimit()
	{
		return ((LocatedParam) _speciesParam).shoveLimit;
	}
	
	/**
	 * \brief Return the shoving interaction distance to be used in shoving
	 * for this species of agent.
	 * 
	 * @return	Double specifying the shoving interaction distance that will
	 * be applied.
	 */
	public Double getInteractDistance()
	{
		return getInteractDistance(this);
	}
	
	/**
	 * \brief Return the shoving interaction distance to be used in shoving
	 * against a specified agent.
	 * 
	 * @return	Double specifying the shoving interaction distance that will
	 * be applied.
	 */
	public Double getInteractDistance(LocatedAgent aLoc)
	{
		return getShoveRadius() + aLoc.getShoveRadius() + getShoveLimit();
	}
	
	/**
	 * \brief Return the fraction of mass that is transferred to the new agent
	 * on cell division.
	 * 
	 * @return	Double stating the fraction of mass that is transferred to the
	 * new agent on cell division.
	 */
	public Double getBabyMassFrac()
	{
		return ExtraMath.deviateFromCV(getSpeciesParam().babyMassFrac,
											getSpeciesParam().babyMassFracCV);
	}
	
	/**
	 * \brief Return the agent radius at which cell division is triggered.
	 * 
	 * @return	Double stating the agent radius at which cell division is
	 * triggered.
	 */
	public Double getDivRadius()
	{
		return ExtraMath.deviateFromCV(getSpeciesParam().divRadius,
											getSpeciesParam().divRadiusCV);
	}
	
	/**
	 * \brief Return the agent radius at which cell death is triggered
	 * 
	 * @return	Double stating the agent radius at which cell death is triggered
	 */
	public Double getDeathRadius()
	{
		return ExtraMath.deviateFromCV(getSpeciesParam().deathRadius,
											getSpeciesParam().deathRadiusCV);
	}
	
	/**
	 * \brief Report if this agent is attached to a surface.
	 * 
	 * @return Boolean noting whether the agent is attached to a surface.
	 */
	public Boolean isAttached()
	{
		return _isAttached;
	}
	
	/**
	 * \brief Return the active fraction of this agent.
	 * 
	 * @return	Double value stating the active fraction of this agent.
	 */
	public Double getActiveFrac()
	{
		return 1.0;
	}
	
	/**
	 * \brief Return the color assigned to this agent in POV-Ray output.
	 * 
	 * TODO Rob 13Mar2015: Consider deleting as part of move away from POV-Ray.
	 * 
	 * @return	Colour assigned to this agent as specified in the protocol
	 * file.
	 */
	public Color getColor()
	{
		return _species.color;
	}

	/**
	 * \brief Return the colour assigned to any capsules contained in this
	 * agent in POV-Ray output.
	 * 
	 * @return	Colour assigned to this agent capsules as specified in the
	 * protocol file.
	 */
	public Color getColorCapsule()
	{
		return Color.green;
	}

	/**
	 * \brief Return the location of this agent.
	 * 
	 * @return	ContinuousVector stating the location of this agent.
	 */
	public ContinuousVector getLocation()
	{
		/*double magnitude = locationHeight.distance(center1);
		if (magnitude == 0)
		{
			randomizeOrientation();
			 magnitude = locationHeight.distance(center1);
		}
		double magX = ((locationHeight.x - center1.x) / magnitude) * (_height) ;
		double magY = ((locationHeight.y - center1.y) / magnitude) * (_height) ;
		double magZ = ((locationHeight.z - center1.z) / magnitude) * (_height) ;
		_location.x=center1.x-magX;
		//_location.y=center1.y-magY;
		_location.z=center1.z-magZ;*/
		return _location;
	}
	
	
	/**
	 * \brief Return the location + height of this  cylindrical agent.
	 * 
	 * @return	ContinuousVector stating the location + height of this cylindrical agent.
	 */
	
	
	
	public ContinuousVector getLocationHeight()
	{
		if(_species.speciesClass.equals("Fungus")) {
//			Double xlocheight = _location.x+ _directionx *_height;
//			/* The 2.0 * Math.random() is inserted to slant the position of the hyphal agent*/
//			/*Double ylocheight =_location.y+2.0*Math.random();*/ 
//			Double ylocheight =_location.y+ (random.nextBoolean() ? 1 : -1 )* _directiony;
//			
			/*Double xlocheight = _location.x+_height;
			Double ylocheight =_location.y;
			locationHeight.set(xlocheight, ylocheight, _location.z);*/
			
			
			
			//To consider the angle and thus torque, the x coordinate is
			//going to be changed as in the idynobacillus code. 
			/*double magnitude = locationHeight.distance(_location);
			if (magnitude == 0)
			{
				randomizeOrientation();
				 magnitude = locationHeight.distance(_location);
			}
			double magX = ((locationHeight.x - _location.x) / magnitude) * (_height) ;
			//double magY = ((locationHeight.y - _location.y) / magnitude) * (_height) ;
			//double magZ = ((locationHeight.z - _location.z) / magnitude) * (_height) ;
			
			Double xlocheight = _location.x+magX;
			Double ylocheight =_location.y;
			//Double xlocheight = _location.x;
			//Double ylocheight =_location.y+magY;
			Double zlocheight = _location.z;
			locationHeight.set(xlocheight, ylocheight, zlocheight);
			//System.out.println(locationHeight);
			center1.x=_location.x-magX;
			//_location.y=center1.y-magY;
			center1.z=_location.z-magZ;*/
			return locationHeight;
		}
		else {
			locationHeight.set(_location); /* If the species is not fungus, it will be spherical */
			return locationHeight;
		}
	}
	
	public ContinuousVector getCenter1()
	{
		if(_species.speciesClass.equals("Fungus")) {
			return center1;
		}
		else {
			center1.set(_location); //For spherical, the center1 is the location itself
			return center1;
		}
	}
	/**
	 * \brief Return the slope of the cylindrical axis assuming its along the z axis
	 * 
 * @return slope
 *  */
	
	/*public Double slopeaxis() {
		Double slope = (locationHeight.y - _location.y)/(locationHeight.x - _location.y);
		return slope;
	}
	*/
	
	
	
	/**
	 * \brief Comparator used by AgentContainer.erodeBorder()
	 * 
	 * Comparator used by AgentContainer.erodeBorder()
	 * @author Rob Clegg
	 */
	public static class detPriorityComparator implements java.util.Comparator<Object>
	{
		@Override
		public int compare(Object b1, Object b2)
		{
			Double f1 = ((LocatedAgent) b1).detPriority;
			Double f2 = ((LocatedAgent) b2).detPriority;
			return (int) Math.signum(f1 - f2);
		}
	}

	/**
	 * \brief Comparator used by AgentContainer.erodeBorder()
	 * 
	 * @author Rob Clegg
	 */
	public static class totalMassComparator implements java.util.Comparator<Object>
	{
		@Override
		public int compare(Object b1, Object b2)
		{
			Double f1 = ((LocatedAgent) b1)._totalMass;
			Double f2 = ((LocatedAgent) b2)._totalMass;
			return (int) Math.signum(f1 - f2);
		}
	}
	
	/**
	 * \brief Return the distance from this agent to a ContinuousVector.
	 * 
	 * @param position	ContinuousVector to find distance to.
	 * @return distance between this agent and cV (assuming cyclic boundaries).
	 */
	public Double getDistance(ContinuousVector position)
	{
		return computeDifferenceVector(position).norm();
	}
	
	/**
	 * \brief Return the distance from this agent to another.
	 * 
	 * @param aLoc	LocatedAgent to find distance to.
	 * @return Distance from this agent to that given (assuming cyclic
	 * boundaries).
	 */
	public Double getDistance(LocatedAgent aLoc)
	{
		if(aLoc.getStringClass().equals("Fungus")) {
			ContinuousVector gd=computeDifferenceAxis(aLoc);
			return gd.norm();
		}else {
		return getDistance(aLoc._location);}
	}

	/**
	 * \brief Set the location of this agent to the supplied continuous vector.
	 * 
	 * @param cc	Location which this agent should be assigned to.
	 */
	public void setLocation(ContinuousVector cc) 
	{
		// In a chemostat set the location of the newborns to zero.
		if ( Simulator.isChemostat )
			_location.reset();
		else
			_location.set(cc);
	}
	
	/**
	 * \brief Set the location of the top of the cylindrical agent to the supplied continuous vector.
	 * 
	 * @param cc	Location of the top of cylindrical agent should be assigned to.
	 */
	public void setLocationHeight(ContinuousVector cc) 
	{
		// In a chemostat set the location of the newborns to zero.
		if ( Simulator.isChemostat )
			locationHeight.reset();
		else
			locationHeight.set(cc);
	}
	/**
	 * \brief Return the continuous vector that states this agents move.
	 * 
	 * @return Continuous vector that states this agents move.
	 */
	public ContinuousVector getMovement()
	{
		return _movement;
	}

	/**
	 * \brief Return the index of the grid on which this agent is placed.
	 * 
	 * @return Integer grid index of where this agent is placed.
	 */
	public int getGridIndex()
	{
		return _agentGridIndex;
	}

	/**
	 * \brief Return the LocatedGroup of agents that are present in the
	 * location where this agent is placed.
	 * 
	 * @return	LocatedGroup containing all agents present in the same grid
	 * space as this agent.
	 */
	public LocatedGroup getGridElement()
	{
		return _agentGrid.getShovingGrid()[_agentGridIndex];
	}
	
	/**
	 * \brief Move this agent to another grid index.
	 * 
	 * @param aGridIndex Grid index in which this agent should now be placed.
	 */
	public void setGridIndex(int aGridIndex)
	{
		_agentGridIndex = aGridIndex;
	}

	/**
	 * \brief Return the domain where this agent is contained.
	 * 
	 * @return The domain where this agent is contained (Domain object).
	 */
	public Domain getDomain()
	{
		return _species.domain;
	}
	//return end point
		public static ContinuousVector RotateVector(double theta, 
				ContinuousVector[] v, ContinuousVector[] orientation)
	    {
			//0 = start, 1 = end of vector
			//normalization of localized euclidean vector
			double norm = orientation[1].distance(orientation[0]);

			double vo_mag_x = (orientation[1].x - orientation[0].x) / norm;
	        double vo_mag_y = (orientation[1].y - orientation[0].y) / norm;
	        double vo_mag_z = (orientation[1].z - orientation[0].z) / norm;

	        double mag_x = v[1].x-v[0].x;
	        double mag_y = v[1].y-v[0].y;
	        double mag_z = v[1].z-v[0].z;
	        Quaternion Q1 = new Quaternion((double)0,mag_x,mag_y,mag_z);
	        
	        Quaternion Q2 = new Quaternion((float)Math.cos(theta / 2),
		            (float)(/*vo.x*/vo_mag_x * Math.sin(theta / 2)),
		            (float)(/*vo.y*/vo_mag_y * Math.sin(theta / 2)),
		            (float)(/*vo.z*/vo_mag_z * Math.sin(theta / 2)));
	       Quaternion conjQ2 = Quaternion.Conjugate(Q2);

	        Quaternion Q3;

	        Q3 = Quaternion.Multiply(Quaternion.Multiply(Q2,Q1),conjQ2);

	        ContinuousVector result = new ContinuousVector(v[0].x + Q3.x, v[0].y + Q3.y, v[0].z + Q3.z);
	        return result;
	    }

}