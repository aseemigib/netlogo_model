# netlogo_model
A NetLogo model of thermal protein aggregation

The NetLogo model of thermal protein aggregation is a simple toy model, developed to study kinetics of thermal aggregation. It is created on the basis of experimental results obtained in thermal aggregation experiments carried out for single protein species. The model extrapolates this single-protein thermal aggregation behaviour to multiple protein species aggregating together as a mixture. The model assumes that proteins of different species interact on a 2D toroidal surface. The total number (population) and types of different protein species that can be modelled in a simulation round have no limit in principle, only in terms of computational power.

The model can be used to-

·        Get a general feel of protein thermal aggregation

·        Check how thermal aggregation behaviour of protein changes by changing values of different parameters in the code.

·        Compare simulation results with experimental results of protein thermal aggregation using various protein systems

·        Generate new hypotheses about thermal aggregation of complex protein mixtures

Assumptions of the model (for parsimony)-

·        Proteins exist in one of two states, either folded or unfolded

·        Only unfolded states can bonded to other proteins and form aggregations, and proteins bonded to other proteins cannot fold

·        Thermal aggregation is a reversible process, as proteins can leave the aggregate

·        All proteins within a species will have equal molecular weight

·        The extent of thermal folding-unfolding of protein is a function of its molecular weight

·        All protein species under inspection will form and lose bonds with the same rates

·        Aggregates will be mixed (many different species can be a part of one single aggregate)

·        Proteins can only form links with physically proximate aggregates. 


Each protein agent is characterized by the following state variables-

·        Folded (green) or unfolded (red)

·        Leader (becomes relevant during aggregation process) [NOTE: Unsure if we should include "leader" as it is a bit of a hack to get the proteins to move together, not really a relevant variable]

·        Species (the relative concentrations are defined on the user interface)

·        Molecular weight (the weight of each protein species is defined in the user interface)

·        DeltaG, Keq, k (thermodynamic parameters that are a function of molecular weight)

The following properties are considered ‘globals’, for being equal for all agents -

·        h (required to obtain parabolic nature in deltaG Vs. temperature graph)

·        scaling factor (required to fix the curvature of parabola in deltaG Vs. temperature graph)

Values for these properties are chosen to match the simulation results to the experimental results as closely as possible. These values could, in theory, be changed in code as per the user’s requirement.

The user specifies the number of proteins in the user interface, along with the number of species and each species' weight and relative abundance. To work, the relative proportions must sum to 1. At initialization, all protein agents begin in the ‘folded’ state, as is typical in experimental settings. The spatial coordinates of the agents are initialized randomly, and the agent size is a function of the species weight. When the simulation is started, the properties of turtles start changing each time-step. The changes happen to positions of turtles and their state (from folded to unfolded). Both motion and folding/unfolding depend on the temperature and molecular weight.

Each protein (or aggregate of proteins) is subject to Brownian motion. Each time-step, the heading is randomly drawn from the uniform distribution, and distance traveled is inverseley proportional to the aggregate weight and linearly proportional to the temperature.

The probability that a protein folds (if unfolded) or unfolds (if folded) is also a function of temperature and molecular weight, as depicted in Figure S3A. The formulae determining these probabilities are coded in terms of deltaG, Keq, and k, defined at initialization as a function of T and molecular weight. 

The reversible process of thermal aggregation will happen based on user-defined rates of aggregation and disaggregation, along with the number of physical interactions. When two (or more) unfolded proteins drift close together, they form a link with a probability specified by the user. Unfolded proteins that form links to other proteins are considered in the aggregated (vs disaggregated) state, which will be decided by rates of aggregation and disaggregation. Each agent that is part of an aggregate has a probability of leaving the aggregate on each time-step, at which point its bonds are destroyed and it breaks free into a soluble phase.

The main order parameters (or "dependent variables") in the model are total percent protein folded, total percent protein soluble, and the distribution of aggregate size. To get these results out of the simulation environment in a processable format, one can use BehaviorSpace, a tool in NetLogo to carry out simulation runs in a high-throughput manner and obtain results in a spreadsheet format. All the simulation results shown in different figures throughout the manuscript are created from the data obtained through BehaviorSpace.

For detailed description of elements in interface, please refer to supplementary information at https://sites.google.com/view/simulation-aseem/home

