# Best practices for computing free energy profiles and potentials of mean force
Best Practices document to be submitted to the Living Journal of Computational Molecular Science

Best practices for computing a potential of mean force and free energy profiles

Contributors: Anthony Hazel, Michael Shirts (mrshirts), Sunny Hwang, Alan Grossfield (agrossfield), Baron Peters (n/a), Heather Mayes (hmayes)
People to invite: Jon Whitmer, Andrew Ferguson, Francesco Gervasio

This is intended to be a simple of workflow for designing a free energy curve / PMF calculation


The Github document is:

https://github.com/team-mayes/best_practice_pmf




#Introduction and Scope
State goal: Help researchers get started productively. The theory is well laid out (will provide references). RC/CV very difficult and system dependent but we can include guidelines.

To accomplish this, we will introduce key terms and some common methods.
Discuss what can do, but that will not be the focus on this document
Rates (need for diffusion coefficient)

We will not cover the basics.  Consult other best practices documents
   - MD basics (TODO: GitHub link)
   - MD setup, biomolecular setup (TODO: GitHub link)
   - Statistical error and uncertainty analysis? (TODO: GitHub link)
 You may want to glance at these first to check if you're comfortable with the material.

We will not cover alchemical methods.  See their best practices document  
   (add GitHub link)


Distinguish between methods that choose a coordinate a priori, and those that do not. We will mention those that do not (MSMs, TPS) but leave outside the scope.
MSMs gives an RC but hard to relate to any theoretical condition, change temp, anything.
What are you going to do: RC for further sampling? If you want to generate low D representation MSMs are great; if want to discover an RC, TPS
In the scope will only be combinations of sampling and constructing PMFs from a priori chosen coordinates. Will discuss a first step of learning about the system and metastable states before choosing:
- US/Replica exchange
- Metadynamics
- adaptive force biasing
- Some on weighted ensemble

There are many other enhanced sampling methods (multicanonical, Wang-Landau, Adaptive Umbrella Sampling, Tsallis statistics).  Can't cover them all.  Suggest reviews?

Discuss analysis of results with committor analysis



#Terminology: Note that terms are, and how they are actually used;
   - Discuss both the definitions and actual use (not always strictly correct!)
   - PMF (analogous to probability density)  vs free energy surface (analogous to probability); Only ABF gives a PMF, umbrella sampling and metadynamics give free energy curves.

PMF is a free energy as a function of almost, but not all the coordinates
 - PMF(x0,x1. . . ) = -kb T log x  x0,x1 exp(- beta U(x)) dx (all but x0, x1)
 - Only defined up to a constant.
 - write equation to show distinction between PMF and free energy
 - Jacobian correction
    - Common jacobians -- distance (J= 4 pi r^2), angles (J = 1)
    - Other common rxn coords J is unknown (native contacts, RMSD)
    - Make a figure showing free energy and PMF for a couple of simple
    systems (2 ideal gas particles, 2 LJ particles?)

PMF is a fundamentally a continuous function (unless the rxn coord is intrinsically discrete), but we have to approximate the continuous function when we estimate it. Two approaches to approximate:
-Expectations of histograms or other indicator functions
Delta A(xi) = \sum_{k=1}^K \sum_{n=1}^N \sum W_i(x_n) I_k(x_n), where I_k is some local indicator function, like a top hat function centered at x_n, or a gaussian centered at x_n.
-Some sort of fit from the data to a continuous function with different parameters, for example by a least squares fit or Kullbeck-Liebler divergence. (refs)  Useful if you're using the curves as input to simpler calculations

What is the reaction coordinate/OP/CV?
   - From the theory of phase transitions:
   - OP distinguishes two metastable states
   - CV any collection of variables
   - RC


# Considerations before you start computing   
Is the question well-posed? Is the project doable?
   - Sturgeon’s law as an organizing principle (“Ninety percent of everything is crud”) -- what can we do to drop the percentage?
   -Recommend starting with conventional MD (possibly with T Rep Ex) to understand the system
   - What is the choice for (x0,x1, . . . ) that you don’t integrate over?
      - Careful choice needed, poor choice uninformative or misleading
      - All choices are true, only some are useful or informative.
         - if you try to interpret mechanism, you're asserting that all other degrees of freedom are fast (and thus in equilibrium) with respect to the chosen ones.
      - What are relevant degree or degrees of freedom?
   - What kinds of barriers are expected along the reaction coordinate?
   - Lengthscale of features on reaction coordinate   
   - Timescales for relaxation of orthogonal degrees of freedom? How long do individual trajectories need to be?  
   - Can we implement the reaction coordinate efficiently as a bias?
      - Most packages implement some simple restraints
      - Plugins: COLVARS, PLUMED
   - How many dimensions?
      - Usually, people only do 1 or 2 dimensions.
      - More dimensions makes it more complicated, both to compute, and more importantly to understand what the heck is going on.  Most of the recommendations will be for 1 and 2 dimensions, with a separate section for ndim > 2.

# Sampling techniques: mention what they do.
   -Restrict ourselves to conventional/T Rep EX, US / HREx , Metadynamics, and adaptive force biasing; add strengths and weaknesses of each method (dynamic range to be surmounted)
   - Conventional MD / T Rep Ex (conventional histogram; MSM to analyze)
   - Direct counting/MSM: Usually bad, because only samples the minima, not the barriers.  Will fail even in the absence of barriers -- just need a significant (> 5 kT, maybe) range of free energies
   - Umbrella sampling / Ham Rep Ex (need weighted)
   - Design simulation so it can be tested
      - Prefer running each window multiple times, construct as independently as possible, best way to get error bars
   - Tests for convergence and self-consistency
      - For Ham Rep Ex 1.	Neale, C., J.C.Y. Hsu, C.M. Yip, and R. Pomès. 2014. Indolicidin binding induces thinning of a lipid bilayer. Biophys. J. 106: L29–31.

# Computing the free energy/PMF from US/Ham Rep Ex
   - WHAM
      - Several implementations: mine (http://membrane.urmc.rochester.edu/content/wham), g_wham in gromacs, built in to CHARMM
      - Special case for extreme range of free energies: https://github.com/dejunlin/gwham
      - Minimal tests you need to do
         - Check that histograms overlap
         - Check that bin size is small enough that it’s not altering the answer
         - Tolerance must be small enough that the PMF isn’t changing
         - Make sure units match
            - Mention Amber issue with angle restraints as example
         - Make sure factor of 1/2 matches
         - This is the terminology from Grossfield's WHAM, look at g_wham, etc and check how they talk about it
      - Computing other properties using umbrella-sampled data
         - In principle, WHAM can do it, but probably not the best choice because it's harder to control histogramming artifacts. Recommend
         MBAR
   - MBAR (https://github.com/choderalab/pymbar), also an implementation in Toni Mey's github repo, add link
      - @mrshirts: insert something smart here. :)
   - TRAM would also be a method to compute a pmf (https://github.com/markovmodel/thermotools/tree/devel/ext) (combine weighted and unweighted)
      - Ask Toni Mey to insert something smart here


- Metadynamics (well-tempered metadynamics)
   - why choose which variant of metadynamics
   - risks of choosing too fact an annealing time in well-tempered
   - Use estimator from bias vs. run for a while, then fix the bias and sampling on that then reweight (in which case you're probably better off just using original metadynamics)
   - One of the people we've discussed inviting should really write this

- Weighted ensemble; good for complex/computationally expensive reaction coordinates
   - explain the principle, link to original Huber and Kim 1996 and Chong and Zuckerman 2017
   - advantage is that you don't need to put the bias into the dynamics
      - great for complex or expensive reaction coordinates
      - gives correct kinetics for free in additional to free energy curve
      - WESTPA as reference implementation
- Adaptive biasing force--samples and computes PMF
   - Extended ABF
   - Do we need an expert on this? Jerome Henin maybe?

#  Analysis common to all methods

- Signatures of what has gone wrong within post-processing/ troubleshooting guide
   - Units (open MM is good, Amber has inconsistent units in angle restraints)
   - Stupid factor of 2 prefactors are inconsistent across codes
   - Error
   - Show some plots of what “bad” data will look like
      - I have some examples on butane I can dig up
- Signatures of what has gone wrong in the simulation itself
   - Hysteresis
   - Slow relaxation -- trim off the front of individual windows, look for systematic change in the resulting curve, cite Neale and Pomes
   - Look for structural similarity in neighboring bins
   - Look for hidden orthogonal barriers in HamRepex by failed replica percolation (Neale and Pomes)
- Committor analysis to verify you've found transition state
- Computing error bars on PMFs / free energy curves
   - Bootstrapping, etc: needs correlation as input, not particularly reliable.  
      - Not necessarily the correlation time of the rxn coord
      - Example: ion through ion channel, z-fluc of ion is faster than protein fluctuations
- Do multiple independent runs
- PMF only defined up to constant, how do you align curves to get error bars
   - Several possible solutions (matrix of differences, bin-to-bin differences) credit Chodera, Louis Smith came up with bin-to-bin diffs, not published




Notes not (yet) incorporated:
Restrict ourselves to a priori choices of RC. Mention that others exist
--approximate continuous function with discrete samples (most of the time)
--Projection onto the state space
--Histograms vs. laying down Gaussians
--Histograms always underestimates barrier heights (not the worst problem)

MSM as an adaptive sampling method (will exclude this from our scope)

An article on “Good ideas that didn’t work” -- examples of problems in PMF calculations
   Suggestion by Dan Zuckerman
   E.g. “I chose the wrong reaction coordinate and it ruined my life”

Weighted ensemble;

Colvers in Plumed

Wang-Landau--MetaD is the continuous version


What are the citations that demonstrate the issues here?
