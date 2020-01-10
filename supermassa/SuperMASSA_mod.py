# Copyright 2012 Oliver Serang
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#
#    Adapted by Gabriel Gesteira (gabrielgesteira@usp.br) in December 16th, 2019
#    to integrate with R/MAPpoly package (https://github.com/mmollina/MAPpoly)

import tempfile
import numpy
import math
import sys
import getopt
import cProfile
import itertools
import copy
import pprint
import matplotlib
import numpy
from matplotlib import pylab as P
# import pylab as P
import sys
import math
import random
import sys

def log(val_f):
    if val_f == 0:
        return -float('inf')
    return math.log(val_f)

def log_pow(base, exponent):
    if exponent == 0:
        return 0
    return exponent*log(base)

def normal_pdf(deviation, sigma):
    return 1/math.sqrt(2*math.pi*sigma)*math.exp(- deviation*deviation / (2*sigma*sigma))

def log_normal_pdf(deviation, sigma):
    return -log(math.sqrt(2*math.pi*sigma)) + -deviation*deviation / (2*sigma*sigma)

#def log_squared_normal_pdf(deviation, sigma):
#    return -2*log(math.sqrt(2*math.pi*sigma)) + -2*deviation*deviation / (2*sigma*sigma)

# fixme: some error here
# def binomial(n,x):
#     if n > x:
#         return 0
#     return prod( [ float(i) for i in xrange(x+1, n+1) ] ) / prod([float(i) for i in xrange(2, n-x+1)])

def log_binomial(n, x):
    if x > n:
        return log(0)
    return sum( [ log(i) for i in xrange(x+1, n+1) ] ) - sum([log(i) for i in xrange(2, n-x+1)])

log_sums_cache = {}
tot = 0.0
for i in xrange(1, 400):
    tot += log(i)
    log_sums_cache[i] = tot
# not accurate for log sums, but for factorial interpretation
log_sums_cache[0] = 0 
#print log_sums_cache

def log_binomial_cached(n, x):
    if n in log_sums_cache:
        if x > n:
            return log(0)
        else:
            # if in the cache, return the value
#            print 'all vals in cache', n, x
#            print log_sums_cache[n] - log_sums_cache[x] - log_sums_cache[n-x]
            return log_sums_cache[n] - log_sums_cache[x] - log_sums_cache[n-x]
    else:
        print 'fixme: implement cache extension'
        print 'wanted log_binomial of', n, x
        raise Exception('fixme: remove this later')
        sys.exit(1)
        # extend the cache
        for i in xrange(max(), max(n,x)+1):
            tot += log(i)
            log_sums_cache[i] = tot

def log_add(logA,logB):
    if logA == log(0):
        return logB
    if logA<logB:
        return log_add(logB,logA)
    return log( 1 + math.exp(logB-logA) ) + logA

def log_sum(lst):
    # given a list of floats [log(a), log(b), log(c), ...]  compute
    # log(a+b+c+d+...)
    
    # note: this could be faster by not calling log_add

    result = -float('inf')
    for i in lst:
        result = log_add(i, result)
    return result

def sample_from_distribution(dist):
    u = random.uniform(0.0, 1.0)
    
    cum = 0.0
    for outcome in dist:
        cum += dist[outcome]
        if cum >= u:
            return outcome
    print 'warning: sum of cdf from sample_from_distribution was < 1.0'

def sample_from_normal(mean, var):
    return random.normalvariate(mean, var)

def log_probability_of_histogram_given_distribution(counts_histogram, probability_dist):
    # this uses a log multinomial
    total_counts = sum( counts_histogram.values() )
    
    result = 0.0
    total_remaining = total_counts
    for g,c in counts_histogram.items():
        if c == 0:
            continue
        result += log_binomial_cached( total_remaining, c )
        if g not in probability_dist:
            return -float('inf')
        result += log_pow( probability_dist[g], c)
        total_remaining -= c
    return result

###### end of numerics.py

def draw_histogram(experimental_geno_hist, theoretical_geno_hist, save_file_prefix = '', file_type = 'png'):
    if len(experimental_geno_hist) == 0:
        raise Exception('cannot draw histogram of no genotypes')

    if save_file_prefix == '':
        P.ion()
        P.figure(0)
    else:
        P.figure(0)
    P.clf()

    total = sum( experimental_geno_hist.values() )

    P.bar( *zip(* [ (g[0], theoretical_geno_hist[g]) for g in sorted(theoretical_geno_hist)] ), width = 1, color = 'blue', alpha = 0.9 , label = 'Expected')
    P.bar( *zip(* [ (g[0], experimental_geno_hist[g] / float(total)) for g in sorted(experimental_geno_hist)] ), width = 1, color = 'red', alpha = 0.7, label = 'Observed' )

    P.xlabel('Doses of allele 1', size='x-large')
    P.ylabel('Frequencies', size='x-large')

    P.legend(loc='upper right')
    
    if save_file_prefix == '':
        P.draw()
    else:
        P.savefig(save_file_prefix + '_hist.' + file_type, dpi=600)

def draw_genotypes(individuals_to_data, individuals_to_genotypes, figure_number, save_file_prefix = '', file_type = 'png'):
    # if set(individuals_to_data) is a superset of
    # set(individuals_to_genotypes), then only draw the points for the
    # assigned genotypes
    if len(individuals_to_genotypes) == 0:
        raise Exception('cannot draw scatterplot of no genotypes')

    base_colors = ['red', 'blue', 'green', 'black', 'pink', 'orange', 'purple', 'yellow', 'cyan', 'magenta', 'gray']
    base_markers = [ 's' , 'o' , '^', '+' , 'd' , 'h' , '>' , 'p' , '<' , 'v' , 'x', '1' , '2' , '3' , '4', 'None' , ' ' , '' , '$...$']

    # make base colors and markers the same length
    base_display_number = min(len(base_colors), len(base_markers))
    base_colors = base_colors[:base_display_number]
    base_markers = base_markers[:base_display_number]

    # pair every color with every marker
    colors = []
    markers = []
    for j in xrange(0, base_display_number):
        for i,c in enumerate(base_colors):
            colors.append(c)
            markers.append(base_markers[(j+i)%base_display_number])

    if save_file_prefix == '':
        P.ion()
        P.figure(figure_number)
    else:
        P.figure(figure_number)
    P.clf()

    if len(set(individuals_to_genotypes.values())) > min(len(colors), len(markers)):
        print 'Error: cannot plot more than', min(len(colors), len(markers)), 'distinct colors at this time (automated color generation could lower contrast between adjacent classes)'
        sys.exit(2)
    
    genos_to_x_and_ys = {}
    for individual, geno in individuals_to_genotypes.items():
        if geno not in genos_to_x_and_ys:            genos_to_x_and_ys[geno] = individuals_to_data[individual]
        else:
            genos_to_x_and_ys[geno].extend(individuals_to_data[individual])

    x_max = max([ max([x_and_y[0] for x_and_y in genos_to_x_and_ys[g]]) for g in genos_to_x_and_ys ])
    y_max = max([ max([x_and_y[1] for x_and_y in genos_to_x_and_ys[g]]) for g in genos_to_x_and_ys ])
    x_or_y_max = max(x_max, y_max)
    for g, c in zip(sorted(genos_to_x_and_ys), colors):
        # make the theoretical series for that genotype
        if g[0] > 0.0:
            r_limit_x = x_max / g[0]
        else:
            r_limit_x = 1.0e10

        if g[1] > 0.0:
            r_limit_y = y_max / g[1]
        else:
            r_limit_y = 1.0e10

        r_limit = min( r_limit_x, r_limit_y )
        s = [ (0,0), (g[0]*r_limit, g[1]*r_limit) ]
        P.plot( *zip(*s),  c = c, linewidth=2)

    for g, c, m in zip(sorted(genos_to_x_and_ys), colors, markers):
        P.scatter( *zip(*genos_to_x_and_ys[g]), c = c , marker = m, s=40)

    # this adds a little buffer for the figure; note that pylab.axes should probably be used somehow
    P.scatter([-0.05*x_max,-0.05*x_max,x_max,x_max],[-0.05*y_max,y_max,-0.05*y_max,y_max], visible = False)
    P.xlabel('Intensity of allele 1', size='x-large')
    P.ylabel('Intensity of allele 2', size='x-large')
    P.xlim(0,x_or_y_max)
    P.ylim(0,x_or_y_max)

    if save_file_prefix == '':
        P.draw()
    else:
        P.savefig(save_file_prefix + '_scatter.' + file_type, dpi=600)

def getchar():
    print ''
    Var = raw_input("press a key")

def genotype_histogram(P, genotypes):
    result = dict( [ ((i, P-i), 0.0) for i in xrange(0, P+1) ] )
    
    for g in genotypes:
        result[g] += 1
    return result

def histogram(dictionary):
    # builds a histogram of the outputs of the dictionary
    result = {}
    for val in dictionary.values():
        if val in result:
            result[val] += 1
        else:
            result[val] = 1
    return result

def nonzero_distribution(dist):
    return dict([(key, value) for key, value in dist.items() if value != 0.0])

def indent(i_val):
    for i in range(0, i_val):
        print '  ',

        ############################## end of utilities.py

# note: performance can be improved by making functions that don't use
# "self" static
class Base:
    def __init__(self, filename, ploidy_range, sigma_range, **Kwargs):
        self.load(filename, **Kwargs)
        self.ploidy_range = ploidy_range
        self.sigma_range = sigma_range

    def load(self, filename, heights_or_areas):
        try:
            # f = open(filename)
            # individuals_to_data = {}
            # all_lines = f.readlines()
            # if len(all_lines) == 0:
            #     raise Exception("cannot load an empty file")
            # for line in all_lines:
            #     line_lst = line.split()
            #     if len(line_lst) == 5:
            #         individual, x_height, y_height, x_area, y_area = line_lst
            #         if heights_or_areas == 'heights':
            #             pair = numpy.array( [float(x_height), float(y_height)] )
            #         elif heights_or_areas == 'areas':
            #             pair = numpy.array( [float(x_area), float(y_area)] )
            #         else:
            #             raise Exception('load_scatterplot needs to specify height or area')
            #     elif len(line_lst) == 3:
            #         if heights_or_areas == 'areas':
            #             raise Exception('heights_or_areas does not affect 3 column data sets')
            #         individual, x, y = line_lst
            #         pair = numpy.array([float(x), float(y)])
            #     else:
            #         raise Exception('cannot load file with number of columns != 3 or 5')

            #     if numpy.linalg.norm(pair) > 0.0:
            #         # add this to the data set
            #         if individual in individuals_to_data:
            #             individuals_to_data[ individual ].append( pair )
            #         else:
            #             individuals_to_data[ individual ] = [ pair ]
            # self.individuals_to_data = individuals_to_data
            self.individuals_to_data = filename
        except Exception as err:
            raise Exception('Error loading file: ' + str(err))        

    def get_mean_for_each_individual(self):
        individuals_to_mean_data = {}
        for individual, data in self.individuals_to_data.items():
            mean = sum( data ) / float(len(data))
            individuals_to_mean_data[individual] = mean
        return individuals_to_mean_data

    def get_individuals_to_log_likelihoods(self, ploidy, sigma, **Kwargs):
        # takes the product of all data points for each genotype state
        # for each individual
        individuals_to_log_likelihoods = {}
        for individual in self.individuals_to_data:
            sum_log_likelihoods = None
            for pair in self.individuals_to_data[individual]:
                log_likelihoods = self.log_likelihoods_of_genotype_states(ploidy, sigma, pair)

                if sum_log_likelihoods == None:
                    sum_log_likelihoods = log_likelihoods
                else:
                    for g in log_likelihoods:
                        sum_log_likelihoods[g] += log_likelihoods[g]
            individuals_to_log_likelihoods[individual] = sum_log_likelihoods
        return individuals_to_log_likelihoods

    def greedy(self, individuals_to_log_likelihoods, allowed_genotypes):
        individuals_to_best_configurations = {}
        overall_sum_log_likelihoods = 0.0
        for individual, log_likelihoods in individuals_to_log_likelihoods.items():
            best_log_likelihood, best_genotype = max([ (log_like, g) for g, log_like in log_likelihoods.items() if g in allowed_genotypes ])
            individuals_to_best_configurations[ individual ] = best_genotype
            overall_sum_log_likelihoods += best_log_likelihood
        return overall_sum_log_likelihoods, individuals_to_best_configurations

    def log_likelihood_of_configuration(self, inds_to_genos, individuals_to_log_likelihoods):
        return sum([ individuals_to_log_likelihoods[ind][geno] for ind, geno in inds_to_genos.items() ])

    def log_likelihood_of_genotype(self, sigma, x_and_y_G, x_and_y):
        ## fixme: this should to include the fact that the observed x
        ## and y values can never be negative
        theoretical_x_and_y_normalized = numpy.array(x_and_y_G) / float(sum(x_and_y_G))
        x_and_y_normalized = x_and_y / float(sum(x_and_y))
        
        dst = numpy.linalg.norm(x_and_y_normalized - theoretical_x_and_y_normalized)
        return log_normal_pdf(dst, sigma )

    def log_likelihoods_of_genotype_states(self, ploidy, sigma, x_and_y):
        dist = {}
        normalization_constant = log(0)
        for x_G in xrange(0, ploidy+1):
            y_G = ploidy - x_G
            geno = (x_G, y_G)
            log_like = self.log_likelihood_of_genotype(sigma, geno, x_and_y)
            dist[ geno ] = log_like
        return dist

    def optimal(self):
        pass
    def optimal_greedy(self):
        pass

    def best_log_score_posterior_configuration_parameters_for_MAPs(self, log_joints_and_inds_to_genotypes_and_parameters):
        # over all parameters and results, for each set of parameters
        # to marginalize, take take the MAP. then marginalize over the
        # parameters to marginalize to get a posterior estimate.

        # build a map of marg vars to results; then the MAP of these
        # results can be taken for each marg var config
        params_to_MAPs = {}
        best_log_joints_and_inds_to_genotypes_and_parameters_for_marg_vars = []
        for log_joint, config, params in log_joints_and_inds_to_genotypes_and_parameters:
            marg_vars_and_outcomes = tuple( [ (param, val) for param, val in params.items() if param not in self.vars_for_MAP ] )
            MAP_vars_and_outcomes = [ (var, params[var]) for var in self.vars_for_MAP ]

            if marg_vars_and_outcomes in params_to_MAPs:
                params_to_MAPs[marg_vars_and_outcomes].append( (log_joint, MAP_vars_and_outcomes, config) )
            else:
                params_to_MAPs[marg_vars_and_outcomes] = [ (log_joint, MAP_vars_and_outcomes, config) ]

        for marg_vars_and_outcomes, results in params_to_MAPs.items():
            log_joint, MAP_vars_and_outcomes, config = max(results)
            params_dict = dict(marg_vars_and_outcomes)
            params_dict.update( dict(MAP_vars_and_outcomes) )
            best_log_joints_and_inds_to_genotypes_and_parameters_for_marg_vars.append( (log_joint, config, params_dict) )

        log_joint_max, inds_to_genotypes_max, parameters_max = max(best_log_joints_and_inds_to_genotypes_and_parameters_for_marg_vars)

        best_posterior = math.exp( log_joint_max - log_sum([log_joint for log_joint, inds_to_genotypes, params in best_log_joints_and_inds_to_genotypes_and_parameters_for_marg_vars]) )
        return log_joint_max, best_posterior, inds_to_genotypes_max, parameters_max

    def naive_high_quality_genotypes(self, naive_posterior_threshold, individuals_to_genotypes, all_parameter_map, prob_geno_fname):
        individuals_to_log_likelihoods = self.get_individuals_to_log_likelihoods(**all_parameter_map)
        # f = open(str(prob_geno_fname), 'w+')  # Commented to avoid creating temp files
        high_qual_inds_to_genos = {}
        flag = 0
        for ind, log_likelihoods in individuals_to_log_likelihoods.items():
            flag = flag + 1
            log_total = log_sum( log_likelihoods.values() )
            log_chosen = log_likelihoods[ individuals_to_genotypes[ind] ]
            if math.exp(log_chosen - log_total) > naive_posterior_threshold:
                high_qual_inds_to_genos[ind] = individuals_to_genotypes[ind]
        #     if flag == 1:
        #        f.write(str(log_likelihoods.keys()))
        #        f.write('\n')
        #     f.write(ind)
        #     f.write( ' ' )
        #     s = str(log_likelihoods.values())
        #     f.write(s)
        #     f.write('\n')
        # f.close()

        # Creating global variable to export likelihoods
        global itll
        itll = individuals_to_log_likelihoods
        
        return high_qual_inds_to_genos

class SharedPloidyInference(Base):
    def __init__(self, filename, ploidy_range, sigma_range, **Kwargs):
        Base.__init__(self, filename, ploidy_range, sigma_range, **Kwargs)
        self.vars_for_MAP = ['sigma']

    def log_parameter_prior(self, ploidy):
        # treats genotypes as uniform given ploidy
        return log(1/float(ploidy+1))

    def optimal(self, display_progress):
        return self.optimal_greedy(display_progress)

    def get_inds_to_genotypes_given_parent_parameters(self, parameter, sigma, **Kwargs):
        geno_set = set(parameter)
        inds_to_log_likes = self.get_individuals_to_log_likelihoods(sigma=sigma, **Kwargs)
        # return the configuration, not the score
        return self.greedy(inds_to_log_likes, geno_set)[1]

    def optimal_greedy(self, display_progress):
        log_joints_and_inds_to_genotypes_and_parameters = []
        for sigma in self.sigma_range:
            for ploidy in self.ploidy_range:
                individuals_to_log_likelihoods = self.get_individuals_to_log_likelihoods(ploidy, sigma)
                log_pr_D_given_G, inds_to_genotypes = self.greedy(individuals_to_log_likelihoods, set([ (i, ploidy-i) for i in xrange(0, ploidy+1) ]) )
                log_joint = log_pr_D_given_G + self.log_parameter_prior(ploidy)
                log_joints_and_inds_to_genotypes_and_parameters.append((log_pr_D_given_G, inds_to_genotypes,  {'sigma':sigma, 'ploidy':ploidy}))
                if display_progress:
                    print sigma, ploidy, log_joint
        return self.best_log_score_posterior_configuration_parameters_for_MAPs(log_joints_and_inds_to_genotypes_and_parameters)

# rewrite functors with two inherited classes from here
class PopulationInference(Base):
    def __init__(self, filename, ploidy_range, sigma_range, vars_for_MAP, **Kwargs):
        Base.__init__(self, filename, ploidy_range, sigma_range, **Kwargs)
        self.vars_for_MAP = set(vars_for_MAP)
    def optimal_greedy(self, display_progress):
        log_joints_and_inds_to_genotypes_and_parameters = []
        for sigma in self.sigma_range:
            for ploidy in self.ploidy_range:
                individuals_to_log_likelihoods = self.get_individuals_to_log_likelihoods(ploidy, sigma)
                for parameter in self.parameter_range_generator(ploidy = ploidy):
                    T = self.T_generator(ploidy, parameter)
                    T_nonzero = nonzero_distribution(T)
                    log_joint, individuals_to_genotypes = self.single_greedy(ploidy, parameter, T_nonzero, individuals_to_log_likelihoods, sigma)
                    log_joints_and_inds_to_genotypes_and_parameters.append((log_joint, individuals_to_genotypes, {'sigma':sigma, 'ploidy':ploidy, 'parameter':parameter}))
                    if display_progress:
                        print sigma, ploidy, parameter, log_joint
                        #log_total = log_likelihoods.values()
                        #print log_total

        return self.best_log_score_posterior_configuration_parameters_for_MAPs(log_joints_and_inds_to_genotypes_and_parameters)

    def single_greedy(self, ploidy, parameter, T_nonzero, individuals_to_log_likelihoods, sigma):
        # get the greedy result for this set of parameters
        greedy_log_pr_D_given_G, greedy_inds_to_genotypes  = self.greedy(individuals_to_log_likelihoods, set(T_nonzero))
        greedy_C = histogram(greedy_inds_to_genotypes)
        greedy_log_pr_C_given_T = log_probability_of_histogram_given_distribution(greedy_C, T_nonzero)
        greedy_log_pr_D_and_C_given_G_and_T = greedy_log_pr_D_given_G + greedy_log_pr_C_given_T
        log_prior = self.log_parameter_prior(ploidy, parameter=parameter, sigma=sigma)
        greedy_log_score = greedy_log_pr_D_and_C_given_G_and_T + log_prior
        return greedy_log_score, greedy_inds_to_genotypes

    def optimal(self, display_progress, epsilon = 0.01):
        if display_progress:
            print 'seeding with greedy search'
        best_log_score, posterior, config, params = self.optimal_greedy(display_progress)
        if display_progress:
            print 'starting optimal search'

        if self.number_of_marginalized_configs > 0:
            threshold_buffer_for_marg = log(epsilon / float(self.number_of_marginalized_configs))
        else:
            # there is only one confuguration to try, so the buffer
            # can be 0
            threshold_buffer_for_marg = 0

        log_joints_and_inds_to_genotypes_and_parameters = []
        for sigma in self.sigma_range:
            for ploidy in self.ploidy_range:
                individuals_to_log_likelihoods = self.get_individuals_to_log_likelihoods(ploidy, sigma)
                for parameter in self.parameter_range_generator(ploidy = ploidy):
                    T = self.T_generator(ploidy, parameter)
                    T_nonzero = nonzero_distribution(T)

                    # note: this unecessarily sorts the individuals
                    # each iteration. it is arranged this way for
                    # simplicity; plus the performance gain would be
                    # small since runtime is dominated by branch and
                    # bound
                    ordered_individuals, ordered_genotypes = self.get_ordered_individuals_and_genotypes(T_nonzero)

                    # get the greedy log joint and configuration
                    greedy_log_joint, greedy_individuals_to_genotypes = self.single_greedy(ploidy, parameter, T_nonzero, individuals_to_log_likelihoods, sigma)

                    # get a different greedy result based on the best
                    # histogram C
                    hist_greedy_log_joint, hist_greedy_individuals_to_genotypes = self.single_histogram_greedy(ploidy, parameter, T_nonzero, individuals_to_log_likelihoods, ordered_individuals, ordered_genotypes, sigma)

                    # use the greedy results to seed the branch and
                    # bound. also, if a known configuration is
                    # sufficiently better, numerically bound.
                    # furthermore, give a small margin of 0.05 to
                    # prevent bounding from numerical error from
                    # preventing any solutions (in case the optimal
                    # isn't much better than the greedy)

                    # since the branch and bound doesn't use prior,
                    # and all of these do use a prior, then subtract
                    # the prior for this set of parameters from the
                    # threshold; it will be the same as adding it into
                    # the branch and bound. we want the set of
                    # parameters with maximum branch and bound result
                    # (including parameter prior)
                    log_prior = self.log_parameter_prior(ploidy, parameter=parameter, sigma=sigma)
                    threshold = max(greedy_log_joint, hist_greedy_log_joint, best_log_score + threshold_buffer_for_marg) - log_prior - 0.05

                    log_joint, individuals_to_genotypes = self.count_branch_and_bound(individuals_to_log_likelihoods, T_nonzero, threshold, ordered_individuals, ordered_genotypes)
                    log_joint += log_prior

                    best_log_score = max(best_log_score, log_joint)

                    if display_progress:
                        print sigma, ploidy, parameter, log_joint

                    log_joints_and_inds_to_genotypes_and_parameters.append((log_joint, individuals_to_genotypes, {'sigma':sigma, 'ploidy':ploidy, 'parameter':parameter}))
        return self.best_log_score_posterior_configuration_parameters_for_MAPs(log_joints_and_inds_to_genotypes_and_parameters)

    def single_histogram_greedy(self, ploidy, parameter, T_nonzero, individuals_to_log_likelihoods, ordered_individuals, ordered_genotypes, sigma):
        hist_greedy_C = self.greedy_best_possible_histogram_for_theoretical(T_nonzero, len(individuals_to_log_likelihoods))
        hist_greedy_inds_to_genotypes = self.get_individuals_to_genotypes_from_histogram(ordered_individuals, ordered_genotypes, hist_greedy_C)
        hist_greedy_log_pr_C_given_T = log_probability_of_histogram_given_distribution(hist_greedy_C, T_nonzero)
        hist_greedy_log_pr_D_given_G = sum([individuals_to_log_likelihoods[ind][geno] for ind, geno in hist_greedy_inds_to_genotypes.items() ])
        log_prior = self.log_parameter_prior(ploidy, parameter = parameter, sigma = sigma)
        hist_greedy_log_score = hist_greedy_log_pr_D_given_G + hist_greedy_log_pr_C_given_T + log_prior
        return hist_greedy_log_score, hist_greedy_inds_to_genotypes

    def greedy_best_possible_histogram_for_theoretical(self, T_nonzero, n):
        C = {}
        # use sorted order so that results aren't stochastic; I don't
        # believe dictionary items are not guaranteed to be in the
        # same order
        for i, g in enumerate(sorted(T_nonzero)):
            if i == len(T_nonzero)-1:
                break
            expected = int(n*T_nonzero[g])
            C[g] = expected
        # set the last one (it's the one that broke the loop)
        C[g] = n - sum(C.values())
        return C

    def get_ordered_individuals_and_genotypes(self, T_nonzero):
        # sort the genotypes lexicographically (after L1 normalization)
        ordered_individuals = [ (tuple(pair/sum(pair)), individual) for individual, pair in self.get_mean_for_each_individual().items() ]
        ordered_individuals = sorted(ordered_individuals)[::-1]
        # sort the theoretical genotypes lexicographically
        ordered_genotypes = sorted( list(T_nonzero) )[::-1]
        return ordered_individuals, ordered_genotypes

    def count_branch_and_bound(self, individuals_to_log_likelihoods, T_nonzero, threshold, ordered_individuals, ordered_genotypes):
        # arrange the log_likelihoods in the same order
        geometrically_sorted_individuals_and_log_likelihoods = [ (individual, individuals_to_log_likelihoods[individual]) for value, individual in ordered_individuals ]

        best_score, genotypes_to_counts = self.count_branch_and_bound_recursive(geometrically_sorted_individuals_and_log_likelihoods, T_nonzero, ordered_genotypes, {}, 0, 0.0, threshold)
        return best_score, self.get_individuals_to_genotypes_from_histogram(ordered_individuals, ordered_genotypes, genotypes_to_counts)
    
    def count_branch_and_bound_recursive(self, geometrically_sorted_individuals_and_log_likelihoods, T_nonzero, ordered_genotypes, genotypes_to_counts, number_assigned, log_mult_to_get_here, threshold):
        n = len(geometrically_sorted_individuals_and_log_likelihoods)
        number_remaining = n - number_assigned
        
        next_genotype_index = len(genotypes_to_counts)
        if next_genotype_index == len(T_nonzero):
            return log_mult_to_get_here, deepcopy(genotypes_to_counts)

        next_genotype = ordered_genotypes[next_genotype_index]
        recursive_scores_and_genotype_counts = []

        for next_number in xrange(0, number_remaining + 1):
            log_mult_for_assignment = self.log_multiplier(geometrically_sorted_individuals_and_log_likelihoods, number_assigned, T_nonzero, next_genotype, next_number, number_remaining)
            genotypes_to_counts[next_genotype] = next_number
            log_upper_bound_to_finish_from_here = self.log_upper_bound_after_assignment(geometrically_sorted_individuals_and_log_likelihoods, genotypes_to_counts, T_nonzero, number_assigned + next_number)

            new_log_mult_to_get_here = log_mult_to_get_here + log_mult_for_assignment
            if new_log_mult_to_get_here + log_upper_bound_to_finish_from_here < threshold:
                continue

            # cannot bound, so recurse
            best_remaining_score_and_genotype_counts = self.count_branch_and_bound_recursive(geometrically_sorted_individuals_and_log_likelihoods, T_nonzero, ordered_genotypes, genotypes_to_counts, number_assigned+next_number, new_log_mult_to_get_here, threshold)
            recursive_scores_and_genotype_counts.append(best_remaining_score_and_genotype_counts)

        # all outcomes of the code require returning, so the genotype
        # added at this level must be removed
        del genotypes_to_counts[next_genotype]

        if len(recursive_scores_and_genotype_counts) == 0:
            # all for this were bounded
            return -float('inf'), {}
        return max(recursive_scores_and_genotype_counts)

    # mathematical tools used by the branch and bound
    def log_multiplier(self, geometrically_sorted_individuals_and_log_likelihoods, number_assigned, T_nonzero, next_genotype, next_number, number_remaining):
        # computes the probability of drawing these individuals from
        # the distribution and producing this data
        log_prob_drawing_from_T = log_binomial_cached(number_remaining, next_number) + log_pow(T_nonzero[next_genotype], next_number)
        relevant_individuals_and_log_likelihoods = geometrically_sorted_individuals_and_log_likelihoods[number_assigned : number_assigned + next_number]
        log_like = sum( [ log_likelihoods[next_genotype] for individual, log_likelihoods in relevant_individuals_and_log_likelihoods ] )
        return log_prob_drawing_from_T + log_like

    def log_upper_bound_after_assignment(self, geometrically_sorted_individuals_and_log_likelihoods, genotypes_to_counts, T_nonzero, number_assigned):
        # compute a tight upper bound on the best remaining path
        log_best_remaining_likelihood = sum([ max( [ log_likelihoods[g] for g in T_nonzero if g not in genotypes_to_counts ] + [ -float('inf') ] ) for individual, log_likelihoods in geometrically_sorted_individuals_and_log_likelihoods[number_assigned:] ])
        # compute the sum of all probabilities not used by the
        # multinomial (take the max with 0.0 in case of very small
        # numerical error)
        remaining_prob = max(1-sum([ T_nonzero[g] for g in genotypes_to_counts]), 0.0)
        log_best_remaining_multinomial_term = log_pow( remaining_prob, len(geometrically_sorted_individuals_and_log_likelihoods) - number_assigned )
        return log_best_remaining_likelihood + log_best_remaining_multinomial_term

    def get_individuals_to_genotypes_from_histogram(self, ordered_individuals, ordered_genotypes, genotypes_to_counts):
        individuals_to_genotypes = {}
        cumulative_assigned = 0
        for genotype in ordered_genotypes:
            if genotype not in genotypes_to_counts:
                continue
            count = genotypes_to_counts[genotype]
            individuals_to_genotypes.update( dict([ (ind, genotype) for pair, ind in ordered_individuals[cumulative_assigned : cumulative_assigned + count] ]) )
            cumulative_assigned += count
        return individuals_to_genotypes

    def parameter_range_generator():
        pass
    def T_generator():
        pass
    def log_parameter_prior():
        pass

class F1Inference(PopulationInference):
    def __init__(self, filename, ploidy_range, sigma_range, **Kwargs):
        PopulationInference.__init__(self, filename, ploidy_range, sigma_range, vars_for_MAP = ['sigma'], **Kwargs)
        self.number_of_marginalized_configs = sum( [ len([self.parameter_range_generator(ploidy = p)]) for p in self.ploidy_range ] ) - 1
    def parameter_range_generator(self, ploidy):
        for p1 in [ (i, ploidy-i) for i in xrange(0, ploidy+1) ]:
            # use symmetry to only check each pair twice by
            # constraining p2 >= p1
            for p2 in [ (i, ploidy-i) for i in xrange(p1[0], ploidy+1) ]:
                yield((p1, p2))
    def log_parameter_prior(self, ploidy, **Kwargs):
        return log(1/float(math.exp(log_binomial_cached(ploidy+2, 2))))
    def T_generator(self, ploidy, parameter):
        g_p1, g_p2 = parameter
        if ploidy % 2 == 1:
            print 'Warning: using hypergeometric distribution of offspring with an odd ploidy'
        # fill the table with every genotype of interest
        result = dict( [ ((i, ploidy-i), 0.0) for i in range(0, ploidy+1) ] )
        # ploidy should be even, so division by 2 is OK
        for x1 in xrange(0, ploidy/2+1):
            log_prob_x1 = log_binomial_cached(g_p1[0], x1) + log_binomial_cached(ploidy - g_p1[0], ploidy/2 - x1) - log_binomial_cached(ploidy, ploidy/2)
            for x2 in xrange(0, ploidy/2+1):
                log_prob_x2 = log_binomial_cached(g_p2[0], x2) + log_binomial_cached(ploidy - g_p2[0], ploidy/2 - x2) - log_binomial_cached(ploidy, ploidy/2)
                prob = math.exp(log_prob_x1 + log_prob_x2)
                offspring_geno = ( x1+x2, ploidy-(x1+x2) )
                result[ offspring_geno ] += prob
        return result

class SlowF1Inference(F1Inference):
    # overrides the branch and bound routine from ParentalInference to
    # go sloooooowwww.
    def __init__(self, filename, ploidy_range, sigma_range, **Kwargs):
        F1Inference.__init__(self, filename, ploidy_range, sigma_range, **Kwargs)

    def count_branch_and_bound(self, individuals_to_log_likelihoods, T_nonzero, threshold, ordered_individuals, ordered_genotypes):
        return self.genotype_branch_and_bound(individuals_to_log_likelihoods, T_nonzero, threshold)

    def genotype_branch_and_bound(self, individuals_to_log_likelihoods, T_nonzero, threshold):
        ind0 = individuals_to_log_likelihoods.keys()[0]

        inds = individuals_to_log_likelihoods.keys()
        log_likes = [ individuals_to_log_likelihoods[i] for i in inds ]

        genos_to_inds = dict([ (g,i) for i,g in enumerate(T_nonzero) ])
        return self.genotype_branch_and_bound_recursive(inds, log_likes, T_nonzero, threshold, genos_to_inds, { tuple([0]*len(T_nonzero)) : (0.0, {}) }, 0)
    
    def genotype_branch_and_bound_recursive(self, index_inds, index_log_likelihoods, T_nonzero, threshold, genos_to_inds, hists_to_log_scores_and_configs, depth):
        if len(hists_to_log_scores_and_configs) == 0:
            return (-float('inf'), {})

        if depth == len(index_inds):
            return max([(log_score + self.log_hist_prob(hist, T_nonzero), config) for hist, (log_score, config) in hists_to_log_scores_and_configs.items()] )

        ind, log_likes = index_inds[depth], index_log_likelihoods[depth]
        max_remaining_after = sum( [ max(log_l.values()) for log_l in index_log_likelihoods[depth+1:] ] )

        # go through hists at this level
        new_hists_to_log_scores_and_configs = {}
        for hist, (log_score, config) in hists_to_log_scores_and_configs.items():
            for g in T_nonzero:
                i = genos_to_inds[g]
                # get new hist
                new_hist = list(hist)
                new_hist[i] += 1
                new_hist = tuple(new_hist)

                # compute the new score of this config
                new_log_score = log_score + log_likes[g]

                if new_log_score + max_remaining_after >= threshold:
                    new_config = deepcopy(config)
                    new_config[ind] = g
                    if new_hist not in new_hists_to_log_scores_and_configs:
                        new_hists_to_log_scores_and_configs[ new_hist ] = (new_log_score, new_config)
                    elif new_hist in new_hists_to_log_scores_and_configs and new_log_score > new_hists_to_log_scores_and_configs[ new_hist ][0]:
                        new_hists_to_log_scores_and_configs[ new_hist ] = (new_log_score, new_config)

        # recurse using this layer
        return self.genotype_branch_and_bound_recursive(index_inds, index_log_likelihoods, T_nonzero, threshold, genos_to_inds, new_hists_to_log_scores_and_configs, depth+1)

    def log_hist_prob(self, hist, T_nonzero):
        genos_to_counts = dict([ (g, hist[i]) for i,g in enumerate(T_nonzero) ])
        return log_probability_of_histogram_given_distribution(genos_to_counts, T_nonzero)

class all_cached_functor:
    def __init__(self, func):
        self.func = func
        self.cache = {}
    def __call__(self, *args):
        if args in self.cache:
            return self.cache[args]
        else:
            ans = self.func(*args)
            self.cache[args] = ans
            return ans

class single_cached_functor:
    def __init__(self, func):
        self.func = func
        self.cache = {}
    def __call__(self, *args):
        if args in self.cache:
            return self.cache[args]
        else:
            if len(self.cache) > 1:
                self.cache = {}
            ans = self.func(*args)
            self.cache[args] = ans
            return ans

class F1AndParentInference(F1Inference):
    def __init__(self, parent_filename, progeny_filename, ploidy_range, sigma_range, **Kwargs):
        F1Inference.__init__(self, progeny_filename, ploidy_range, sigma_range, **Kwargs)
        self.parent_inference = SharedPloidyInference(parent_filename, ploidy_range, sigma_range, **Kwargs)
        self.number_of_marginalized_configs = sum( [ len([self.parameter_range_generator(ploidy = p)]) for p in self.ploidy_range ] ) - 1
        # note: using this type of cache may increase memory usage significantly
        self.cached_parent_inds_to_log_likelihoods = single_cached_functor(self.parent_inference.get_individuals_to_log_likelihoods)
    def parameter_range_generator(self, ploidy):
        # since the parents are labeled (unlike in the F1Inference
        # model), all (ploidy+1)*(ploidy+1) configurations must be
        # tried
        for p1 in [ (i, ploidy-i) for i in xrange(0, ploidy+1) ]:
            for p2 in [ (i, ploidy-i) for i in xrange(0, ploidy+1) ]:
                yield((p1, p2))
    def log_parameter_prior(self, ploidy, parameter, sigma):
        # sort the parent names just in case they come out in a
        # different order sometimes
        parents_to_genos = dict(zip(sorted(self.parent_inference.individuals_to_data.keys()), list(parameter)))
        # note: this would be much more efficient to use a cached functor for the parent log likelihoods
        #log_likelihood_parents = self.parent_inference.log_likelihood_of_configuration(parents_to_genos, self.parent_inference.get_individuals_to_log_likelihoods(ploidy, sigma))
        like_table = self.cached_parent_inds_to_log_likelihoods(ploidy, sigma)
        log_likelihood_parents = self.parent_inference.log_likelihood_of_configuration(parents_to_genos, like_table)
        # note: calling F1Inference log_parameter_prior is not
        # appropriate because it should return 1/(n+1)(n+1)
        #return F1Inference.log_parameter_prior(self, ploidy) + log_likelihood_parents
        return 2*log(ploidy+1) + log_likelihood_parents

class HWInference(PopulationInference):
    def __init__(self, filename, ploidy_range, sigma_range, **Kwargs):
        PopulationInference.__init__(self, filename, ploidy_range, sigma_range, vars_for_MAP = ['sigma', 'parameter'], **Kwargs)
        self.number_of_marginalized_configs = len(ploidy_range) - 1
    def parameter_range_generator(self, ploidy):
        # note: this doesn't use ploidy. it's accepted as an argument
        # because it needs to match the
        # unique_parents_for_ploidy_generator function.

        return numpy.arange(0.01, 0.99, 0.05)
    def T_generator(self, ploidy, parameter):
        x_freq = parameter
        # fill the table with every genotype of interest
        result = dict( [ ((i, ploidy-i), 0.0) for i in range(0, ploidy+1) ] )
        for g in result:
            x, y = g
            log_prob = log_binomial_cached(ploidy, x) + log_pow(x_freq, x) + log_pow(1-x_freq, ploidy-x)
            result[g] = math.exp(log_prob)
        return result
    def log_parameter_prior(self, ploidy, **Kwargs):
        return log(1 / float(ploidy+1))

def split_range_arguments(str_range, tuple_of_allowed_sizes):
    # accepts things of the form low:high:step, low:high, and val
    str_options = str_range.split(':')
    if len(str_options) not in tuple_of_allowed_sizes:
        raise Exception('this range must specify number of arguments in the following' + str(tuple_of_allowed_sizes))

    if len(str_options) == 3:
        low, high, step_size = str_options
    elif len(str_options) == 2:
        low, high = str_options
        step_size = 2
    elif len(str_options) == 1:
        low = high = str_options[0]
        step_size = 1
    # the values should be numerical, so float will handle any type of
    # result
    low, high, step_size = float(low), float(high), float(step_size)
    if low > high:
        raise Exception('Error: range requires low < high')
    return low, high, step_size

# main and profiling functions
def real_main(options_map):     # Changed argv to options_map
    # try:
    #     opts, args = getopt.getopt(argv, '', ['inference=', 'optimal_inference', 'file=', 'ploidy_range=', 'sigma_range=', 'draw_genotypes', 'print_genotypes', 'print_genotypes_histogram', 'draw_genotype_histogram', 'heights_or_areas=', 'naive_posterior_reporting_threshold=', 'display_progress', 'save_figures=', 'f1_parent_data=', 'slow_inference', 'save_geno_prob_dist='])
    # except getopt.GetoptError as error:
    #     print 'Command line error:' + str(error)
    #     sys.exit(2)

    # options_map = dict( [ (flag[2:], value) for flag, value in opts ] )

    if 'file' not in options_map:
        print 'Error: you must specify a \"file\"'
        sys.exit(2)
    fname = options_map['file']

    if 'save_geno_prob_dist' not in options_map:
        print 'Error: you must specify a \"file\" to save the probability distribution of the genotypes'
        sys.exit(2)
    prob_geno_fname = options_map['save_geno_prob_dist']
        
    if 'inference' not in options_map:
        print 'Error: must specify a type of inference'
        sys.exit(2)
        
    else:
        inference = options_map['inference']
        valid_inference_set = ('ploidy', 'f1', 'hw')
        if inference not in valid_inference_set:
            print 'Error: inference must be one of these', valid_inference_set
            sys.exit(2)

    if 'ploidy_range' in options_map:
        str_ploidy_range = options_map['ploidy_range']
        low_ploidy, high_ploidy, step_size = split_range_arguments(str_ploidy_range, (1,2,3))
        low_ploidy, high_ploidy, step_size = int(low_ploidy), int(high_ploidy), int(step_size)
        ploidy_range = range( low_ploidy, high_ploidy+1, step_size )
        if len(ploidy_range) == 0:
            print 'Error: ploidy range specified is empty'
            sys.exit(2)
    else:
        ploidy_range = range(2, 16+1, 2)

    if 'sigma_range' in options_map:
        str_sigma_range = options_map['sigma_range']
        low_sigma, high_sigma, step_size = split_range_arguments(str_sigma_range, (1,3))
        low_sigma, high_sigma, step_size = float(low_sigma), float(high_sigma), float(step_size)
        # add a small constant in case the low and high are equal
        sigma_range = P.arange( low_sigma, high_sigma+1e-5, step_size )
        if len(sigma_range) == 0:
            print 'Error: sigma range specified is empty'
            sys.exit(2)
    else:
        sigma_range = [ 0.01, 0.02, 0.04, 0.08, 0.16, 0.32]

    if 'heights_or_areas' in options_map:
        if options_map['heights_or_areas'] not in ('heights','areas'):
            print 'Error: heights_or_areas must be either "heights" or "areas"'
            sys.exit(2)
        heights_or_areas = options_map['heights_or_areas']
    else:
        heights_or_areas = 'heights'

    if 'save_figures' in options_map:
        if 'draw_genotypes' not in options_map and 'draw_genotype_histogram' not in options_map:
            print 'Error: cannot save figures if there is nothing to draw'
            sys.exit(2)
        figure_prefix = options_map['save_figures']

    if 'display_progress' in options_map:
        display_progress = True
    else:
        display_progress = False

    if 'f1_parent_data' in options_map:
        if inference != 'f1':
            print 'Error: you must run f1 inference to use f1 parent data'
            sys.exit(1)

    if inference == 'f1':
        if 'f1_parent_data' in options_map:
            parent_fname = options_map['f1_parent_data']
            infer = F1AndParentInference(parent_fname, fname, ploidy_range, sigma_range, heights_or_areas = heights_or_areas)
        else:
            if 'slow_inference' in options_map:
                infer = SlowF1Inference(fname, ploidy_range, sigma_range, heights_or_areas = heights_or_areas)
            else:
                infer = F1Inference(fname, ploidy_range, sigma_range, heights_or_areas = heights_or_areas)
        if 'optimal_inference' in options_map:
            log_joint, posterior, inds_to_genos, parameters = infer.optimal(display_progress)
        else:
            log_joint, posterior, inds_to_genos, parameters = infer.optimal_greedy(display_progress)

        if 'naive_posterior_reporting_threshold' in options_map:
            naive_reporting_thresh = float(options_map['naive_posterior_reporting_threshold'])
            inds_to_genos = infer.naive_high_quality_genotypes(naive_reporting_thresh, inds_to_genos, parameters, prob_geno_fname)

        ploidy, sigma, parents = parameters['ploidy'], parameters['sigma'], parameters['parameter']
        if 'draw_genotype_histogram' in options_map:
            hypergeom = infer.T_generator(ploidy, parents)
            observed_histogram = histogram(inds_to_genos)
            if 'save_figures' in options_map:
                draw_histogram(observed_histogram, hypergeom, figure_prefix)
            else:
                draw_histogram(observed_histogram, hypergeom)

    elif inference == 'hw':
        infer = HWInference(fname, ploidy_range, sigma_range, heights_or_areas = heights_or_areas)
        if 'optimal_inference' in options_map:
            log_joint, posterior, inds_to_genos, parameters = infer.optimal(display_progress)
        else:
            log_joint, posterior, inds_to_genos, parameters = infer.optimal_greedy(display_progress)

        if 'naive_posterior_reporting_threshold' in options_map:
            naive_reporting_thresh = float(options_map['naive_posterior_reporting_threshold'])
            inds_to_genos = infer.naive_high_quality_genotypes(naive_reporting_thresh, inds_to_genos, parameters, prob_geno_fname)

        ploidy, sigma, allele_freq = parameters['ploidy'], parameters['sigma'], parameters['parameter']
        if 'draw_genotype_histogram' in options_map:
            hw_dist = infer.T_generator(ploidy, allele_freq)
            observed_histogram = histogram(inds_to_genos)
            if 'save_figures' in options_map:
                draw_histogram(observed_histogram, hw_dist, figure_prefix)
            else:
                draw_histogram(observed_histogram, hw_dist)

    elif inference == 'ploidy':
        if 'draw_genotype_histogram' in options_map:
            print 'Error: genotype histogram is not relevant w/o parental inference'
            sys.exit(2)

        infer = SharedPloidyInference(fname, ploidy_range, sigma_range, heights_or_areas = heights_or_areas)
        log_joint, posterior, inds_to_genos, parameters = infer.optimal(display_progress)

        if 'naive_posterior_reporting_threshold' in options_map:
            naive_reporting_thresh = float(options_map['naive_posterior_reporting_threshold'])
            inds_to_genos = infer.naive_high_quality_genotypes(naive_reporting_thresh, inds_to_genos, parameters, prob_geno_fname)

# Commented printed outs
# # fixme: only for CGI
#     print 'Results:'
#     print 'log Pr(D, G = g*, parameter*) =', log_joint, ', '
#     print 'Pr(G = g*, parameter* | D) =', posterior, ', '
#     print 'parameter*', parameters
# # fixme: only for CGI
#    # print ''
#     print '\n'
#     if 'print_genotypes' in options_map:
#         # fixme: only for CGI
#         # print '\n'
#         print 'Predicted genotypes:'
#         for i,g in sorted(inds_to_genos.items()):
# # fixme: only for CGI
# #            print i,' ',g
#             print i,':',g," "
#             # fixme: only for CGI
#     # print '\n'
#     if 'print_genotypes_histogram' in options_map:
#         print histogram(inds_to_genos)
#     if 'draw_genotypes' in options_map:
#         if options_map['inference'] == 'f1' and 'f1_parent_data' in options_map:
#             parent_filename = options_map['f1_parent_data']
#             par_infer = SharedPloidyInference(parent_filename, ploidy_range, sigma_range, heights_or_areas = heights_or_areas)
#             par_inds_to_genos = par_infer.get_inds_to_genotypes_given_parent_parameters(**parameters)
#             if 'save_figures' in options_map:
#                 draw_genotypes(par_infer.individuals_to_data, par_inds_to_genos, 2, figure_prefix + '_parents')
#             else:
#                 draw_genotypes(par_infer.individuals_to_data, par_inds_to_genos, 2)
            
#         if 'save_figures' in options_map:
#             draw_genotypes(infer.individuals_to_data, inds_to_genos, 1, figure_prefix)
#         else:
#             draw_genotypes(infer.individuals_to_data, inds_to_genos, 1)

#     if 'save_figures' not in options_map and ( 'draw_genotypes' in options_map or 'draw_genotype_histogram' in options_map ):
#         getchar()

    # Returning objects
    return {'genos': inds_to_genos, 'loglik': itll,  'log.joint': log_joint, 'posterior': posterior, 'args': parameters}

def profile_main(options_map):  # Changed argv to options_map
    cProfile.run('real_main(' + str(options_map) + ')', sort='cumulative')

def folder_prefix(s):
    return s[s.rfind('/', 0, s.rfind('/')-1)+1:]

    class cgi_tools:
        def __init__(self, reform = False):
            self.reform = reform
        def myprint(self, *args):
            if not self.reform:
                print(args)
            for i,x in enumerate(args):
                if type(x) == str:
                    self.myprint(x.replace('\n', '<BR>')),
                elif True: #x is an atom (i.e. int, float, etc.):
                    print x,
                else:
                    self.myprint(x),
                if i != len(args)-1:
                    print ' ',
            print ''
            
def cgi_main(argv):
    try:
        if '--draw_genotypes' in set(argv) or '--draw_genotype_histogram' in set(argv):
            temp_prefix='/gt4sp_1/gdesiqu/'
            temp, temp_fname = tempfile.mkstemp(dir=temp_prefix)
            argv.extend(['--save_figures', temp_fname])
            local_temp_fname = folder_prefix(temp_fname)
        real_main(argv)
        print '<table>'
        if '--draw_genotypes' in set(argv):
            print '<tr><td><h1>Population data</h1><td>'
            print '<img src="/' + local_temp_fname + '_scatter.png" width=600>'
            if '--f1_parent_data' in set(argv):
                print '<tr><td><h1>Parent data</h1><td>'
                print '<img src="/' + local_temp_fname + '_parents_scatter.png" width=600>'
        if '--draw_genotype_histogram' in set(argv):
            print '<tr><td><h1>Genotype frequencies</h1><td>'
            print '<img src="/' + local_temp_fname + '_hist.png" width=600>'
        print '</table>'
    except Exception as error:
        print 'CGI error: ' + str(error)
        sys.exit(2)

# if __name__ == '__main__':
# #    profile_main(sys.argv[1:])
# #    real_main(sys.argv[1:])
#     real_main(sys.argv[1:])
