![alt text](https://github.com/sechseck/epidemics/blob/master/SIR_Figure_1.png "SIR graphs")

# Estimating the SIR epidemic model from observed case numbers

Calculates the expected total case numbers, number of infected over time and number of initially susceptible individuals. It estimates the parameters by solving the SIR differential equations and then using an optimization procedure to find the set of parameters that fit the data observations best.

This approach has been used to model the 2020 cov-19 outbreak in China to good effect up to March and should be applicable to other countries as well. Check the reference section for more information. Everything is provided "as is" and comes without any warranty, to the extent permitted by applicable law. Use at your own risk.
 

## Getting Started
Instructions to get the project up and running on your local machine for data analysis and development.

### Prerequisites
* python 3.7 tested. other versions might work as well
* scipy >= 1.4.1 important, will fail if lower versions are used
* pandas
* numpy
* matplotlib.pyplot


### Installing

put the source file and the case_numbers.csv in the same directory, e.g.
```
$ls
corona_sir.py
case_numbers.csv
```

## Running the test

run corona_sir.py using python and the original case_numbers.csv.
You should get graph and console output.

```
python corona_sir.py
```
The numbers from the console should be equal or very similar to:
```
{   'I0': 45,
    'R0': 0,
    'S': 510528.35167517845,
    'beta': 2.5199541786810506,
    'gamma': 2.3057306408837537}
```
If you get these results, you have a working setup. Congratulations!

## Interpretation of the results
* I0: number of infected individuals at day 0, taken from case_numbers.csv
* R0: number of recovered individuals at day 0 (zero for all intents and purposes)
* S: estimated number of susceptible individuals, i.e. able to become infected due to proximity and susceptibility
* beta: transition rate S->I
* gamma: transition rate I->R

remark: the model does not distinguish between successful recovery or a fatality. 'recovered' is just any live or dead individual that is immune to infection. A first indication is therefore the number of susceptible individuals at the end of the time horizon, as these are guaranteed to be alive.

### Graphs
The unit of the x-Axes of the graphs is number of days as given in case_numbers.csv.

The unit of the y-Axes is the number of individuals in either S,I or R state.


## Extending the model to new countries and updating data

### Adding more case numbers

The first line of case_numbers.csv needs to contain the line
```
country;date;cases
```
add data at the end using:
```
china;17.01.2020;45
china;18.01.2020;62
china;19.01.2020;121
```

### Adding new regions
Similar to adding case numbers, add lines with the new region to case_numbers.csv.
```
new_region;17.01.2020;45
new_region;18.01.2020;62
new_region;19.01.2020;121
```
In corona_sir.py, change 
'''
country = 'china'
'''
to
'''
country = 'new_region'
'''
to the region you just introduced to analyze the new region


## Code available at

The latest version is available from [github](https://github.com/sechseck/epidemics). 

## Authors

**sechseck:** [sechseck](https://github.com/sechseck)

## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/sechseck/epidemics/blob/master/license.md) file for details

## References

**Batista, Milan** (2020). Estimation of the final size of the coronavirus epidemic by the SIR model. [link](https://www.researchgate.net/publication/339311383_Estimation_of_the_final_size_of_the_coronavirus_epidemic_by_the_SIR_model)

https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model

## Acknowledgments

* A big Thank You goes to Milan Batista from Ljubljana for his excellent article. Thanks a lot, hvala lepa!
