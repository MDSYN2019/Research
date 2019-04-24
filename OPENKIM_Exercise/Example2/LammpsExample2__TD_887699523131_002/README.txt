Author: Daniel S. Karls (karl0100 |AT| umn DOT edu), University of Minnesota
Date: 10/10/2014
Version: 002

Revision Notes (version 001 -> 002):
1. Wyckoff information is now reported in the results
2. The use of `sed -i` has been changed to `sed -e` to increase portability
_______________________________________________________________________________________

*                       Section 1: Summary of input variables                         *
_______________________________________________________________________________________

In this Test Driver (TD), the user specifies the following variables through stdin using a Jinja template file named 'pipeline.stdin.tpl' contained in the directory of a Test which uses this TD (the sequence in which they are specified must be the same as below, with 'modelname' being the first):

- 'modelname'
Type: string
Description: Contains the extended-id of a valid KIM Model.  This can be specified by '@< MODELNAME >@' in pipeline.stdin.tpl.  When the pipeline mates a Model with the Test, it will automatically replace '@< MODELNAME >@' with the Model's extended-ID.

- 'element'
Type: string
Description: Contains the abbreviation for the element from which the lattice will be constructed (e.g. Si, Au, Al, etc.).

- 'mass'
Type: double
Description: The molar mass (in g/mol) of the element specified by 'element', e.g. 39.948 for Argon.

- 'latticetypeinput'
Type: string
Description: The type of cubic lattice to construct.  This can be any of "bcc", "fcc", "sc", or "diamond".  Note that this is completely case-insensitive, e.g. "FcC" and "fCC" are both valid specifiers.

- 'a_0'
Type: double
Description: The lattice constant about which most of the lattice spacings will be sampled.  In most cases, this value will be approximately equal to the equilibrium lattice constant and the user will want many samples around this lattice spacing in order to compute the bulk modulus.  Note that the cohesive energy will always be computed at a_0, regardless of what values are used for 'a_min_frac', 'a_max_frac', 'N_lower', 'N_upper', 'samplespacing_lower', and 'samplespacing_upper'.  From a_0, we compute a_min=a_min_frac*a_0 and a_max=a_max_frac*a_0, which will be the minimum and maximum lattice spacings which will be sampled, respectively.

- 'a_min_frac'
Type: double
Description: A fraction which indicates the smallest lattice constant for which the Test will attempt to compute the cohesive energy.  The smallest lattice spacing will be calculated as a_min = a_min_frac*a_0.  Must be strictly greater than zero and strictly less than one.

- 'a_max_frac'
Type: double
Description: A fraction which indicates the largest lattice constant for which the Test will attempt to compute the cohesive energy.  The largest lattice spacing will be calculated as a_max = a_max_frac*a_0.  Must be strictly greater than one.

- 'N_lower'
Type: integer
Description: The number of lattice spacings sampled which are in the interval [a_min, a_0) (where a_min = a_min_frac*a_0).

- 'N_upper'
Type: integer
Description: The number of lattice spacings sampled which are in the interval (a_0, a_max] (where a_max = a_max_frac*a_0).

- 'samplespacing_lower'
Type: double
This parameter controls the density of lattice spacings sampled between a_min and a_0 and must be strictly greater than 1.  Higher values tend to concentrate more sample points near the lattice constant a_0, whereas values near 1 result in a near-uniform sampling of the interval [a_min,a_0].  A value of approximately 5 should be reasonable for most applications, depending on N_lower, a_min, and a_0.

- 'samplespacing_upper'
Type: double
This parameter controls the density of lattice spacings sampled between a_0 and a_max and must be strictly greater than 1.  Higher values tend to concentrate more sample points near the lattice constant a_0, whereas values near 1 result in a near-uniform sampling of the interval [a_0,a_max].  A value of approximately 50 should be reasonable for most applications, depending on N_upper, a_0, and a_max.

_______________________________________________________________________________________

* Section 2: Regarding the meaning of 'samplespacing_lower' and 'samplespacing_upper' *
_______________________________________________________________________________________

Since the majority of lattice spacings sampled are desired to be close to a_0 rather than a_min or a_max, a logarithmic scale is used.
Consider the interval [a_0, a_max] and suppose you have two positive doubles, alpha and beta (beta > alpha). It is possible to construct a logarithmic scale, x, on [0,1] by first breaking up the interval into N equally sized segments of width w through

    w = (beta-alpha)/N

and performing the following loop:

    for(i=0; i<=N; ++i){
        x(i) = ( log(alpha+i*w) - log(alpha) )/( log(beta) - log(alpha) )
    }

However, substituting w and making use of some log identities, we can write this as

    for(i=0; i<=N; ++i){
        x(i) = log( 1+ i*((beta/alpha)/N - i/N )/log(beta/alpha)
    }

Thus, all that matters is the ratio beta/alpha.  Therefore, we can choose to simply fix alpha=1 and specify beta > 1.  The parameter 'beta' is precisely what 'samplespacing_upper' is.  This logarithmic scale on [0,1] is then rescaled to go from a_0 to a_max, giving the following loop to compute the lattice spacings which are sampled in [a_0, a_max]:

    w_upper = (samplespacing_upper - 1)/N_upper

    for(i=0; i<=N_upper; ++i){
        x_upper(i) = a_0 + (a_max-a_0)*[ log( 1 + i*(samplespacing_upper-1)/N_upper )/log(beta) ]
    }

A similar procedure is followed to generate the lattice spacings which are sampled in [a_min, a_0], only with an independent scaling parameter samplespacing_lower.

_______________________________________________________________________________________

*                          Section 3: Invocation of LAMMPS                            *
_______________________________________________________________________________________

Using these input parameters, the TD generates an array of lattice spacings at which the cohesive energy will be computed using the procedure described in Section 2.  Note that there will be a total of numberspacings:=(N_lower + 1 + N_upper) lattice spacings sampled.  N_lower of these spacings will be in [a_min, a_0), one of which is a_min.  N_upper of these spacings will be in (a_0,a_max], one of which is a_max.  Finally, one additional lattice spacing will correspond to a_0.

Next, `sed` is used to replace placholder strings in the 'lammps.in.template' input file script located in the TD directory.  These strings contain
- The atomic symbol
- The atomic mass
- The type of cubic lattice
- The number of lattice spacings that will be sampled
- The particular lattice spacings which will be sampled

The resulting LAMMPS input file is written to 'lammps.in' in the Test Result directory (referred to as "output/" in the TD source).  The LAMMPS script itself builds a single conventional (orthogonal) unit cell with periodic boundary conditions using the LAMMPS 'lattice' command and the 'latticetypeinput' specified by the user.  Notice that a loop over the lattice spacings is natively implemented within lammps.in using the 'next' command (which cycles over the lattice constants specified) and the 'jump' command (which reloads the input script with the new lattice constant).  Running the script creates a single log file, directed to 'lammps.log' in the Test Result directory, which contains the cohesive energy for each of the lattice spacings specified.  This file is then parsed with `grep` and an array containing the cohesive energies corresponding to each lattice spacing is constructed.

Finally, the 'results.edn.tpl' template file is filled in using `sed` to render a valid EDN file named 'results.edn' in the Test Result directory.

_______________________________________________________________________________________

*      Section 4: Using the `testgenie` utility to generate Tests from this TD        *
_______________________________________________________________________________________

This TD contains a collection of template files in a subdirectory named "test_template" which can be used with the `testgenie` KIM utility to generate Tests which use this TD.  The file test_generator.json specifies which Tests to construct by creating JSON dictionaries.  For example, the first JSON dictionary (line 1 in test_generator.json) specifies a Test for fcc Argon which will have a KIM short-id of "TE_778998786610" and uses a_0=5.28, a_min=4.22444, a_max=10.56, N_lower=10, N_upper=18, samplspacing_lower=3, and samplespacing_upper=500.  Using this JSON dictionary, `testgenie` uses the template files in the test_template directory (which are written as Jinja templates) and replaces instances of the specified keywords with these values, after which it creates a directory for the Test in a destination which is specified as an argument to `testgenie`.  For example, to create the two Tests specified by test_generator.json in the directory of this TD, one can issue the following command from the TD root directory:

testgenie --destination ~/openkim-repository/te/ LammpsExample2__TD_887699523131_002

Enter `testgenie --h` to receive more information about `testgenie`.
