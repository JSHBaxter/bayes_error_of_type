# bayes_error_of_type
Exact Bayesian model for errors-of-type with an unknown number of distractors

A model for errors-of-type given an unknown number of distractors is given
in the bayes.py file using the annotators class. To create this model, create
an instance of this class using the number of annotators which will remain 
fixed. The model will then create the cases (available as the cases attribute)
and the basic probability distribution for each case. Observations can be added
one at a time or multiplicitously using the add_case() function. 

To get the results of the model, using the prob() function. It requires taking
a specific regularisation parameter, z, and some description of p and n:
- p can be either: None, (i.e. marginalise over p to get P(n|z,[O]));
                   a rational between 0 and 1 (i.e. evaluate for this p);
                   or "all", (i.e. don't evaluate or marginalise).
- n can be either: None, (i.e. marginalise over n to get P(p|z,[O]));
                   or a positive integer (i.e. evaluate for this n).

prob_mean_p(z,power) gets the expected value of p^power in order to get things
like the mean and variance of p. Lastly, get_mode_p() finds the Gaussian mixture
that approximately fits the distribution P(p|z,[O]).

The poly.py, betalist.py, and rational.py are the exact representation of the
model. The first is a basic class implementing polynomials of multiple
variables with coefficients that are either rational numbers or lists of beta
distributions, hence the remaining two files.

For an example of how to use the model, see the simulation.py file or the
process_tms_data.py file.
