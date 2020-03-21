![alt text](https://github.com/sechseck/epidemics/blob/master/SIR_Figure_1.png "SIR graphs")

# Estimating the SIR epidemic model from observed case numbers Corona / Covid-19
![alt text](https://github.com/sechseck/epidemics/blob/master/PvsO_Figure_2.png "goodness of fit")
<br/><br/>
Calculates the expected total case numbers, number of infected over time and number of initially susceptible individuals. It estimates the parameters by solving the SIR differential equations and then using an optimization procedure to find the set of parameters that fit the data observations best.

This approach has been used to model the 2020 cov-19 outbreak in China to good effect up to March and should be applicable to other countries as well. Check the reference section for more information. Everything is provided "as is" and comes without any warranty, to the extent permitted by applicable law. Use at your own risk.


## Getting Started
Instructions to get the project up and running on your local machine for data analysis and development.

### Prerequisites
* python 3.7 tested. other versions might work as well
* scipy >= 1.4.1 important, will fail if lower versions are used
* pandas 0.25.3 tested
* numpy 1.18.1 tested
* matplotlib.pyplot  3.2.0rc1 tested


### Installing

Put the source file and the case_numbers.csv in the same directory, e.g.
```
$ls
corona_sir.py
case_numbers.csv
```

## Running the test

Run corona_sir.py using python and the original case_numbers.csv.
You should get graph and console output.

```
python corona_sir.py
```
The numbers from the console should be equal or very similar to:
```
{   'I0': 86,
    'R0': 0,
    'S': 261177.81657711515,
    'beta': 1.388808907936211,
    'gamma': 1.1603442769143335}
```
If you get these results, you have a working setup. Congratulations!

## Interpretation of the results
* I0: number of infected individuals at day 1, taken from case_numbers.csv
* R0: number of recovered individuals at day 1 (zero for all intents and purposes)
* S: estimated number of susceptible individuals, i.e. able to become infected due to proximity and susceptibility
* beta: transition rate S->I
* gamma: transition rate I->R

Remark: The model does not distinguish between successful recovery or a fatality. 'recovered' is just any live or dead individual that is immune to infection. A first indication is therefore the number of susceptible individuals at the end of the time horizon, as these are guaranteed to be alive.

### Graphs
The unit of the x-Axes of the graphs is number of days as given in case_numbers.csv.

The unit of the y-Axes is the number of individuals in either S,I or R state.


## Extending the model to new countries and updating data

### Adding more case numbers

The first line of case_numbers.csv needs to contain the line
```
region;date;cases
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
The region to analyze currently defaults to 'china'. This can be changed by using the -r commmand line argument. So to select the case numbers for 'new_region', run
```
python corona_sir.py -r new_region
```

## Code available at

The latest version is available from [github](https://github.com/sechseck/epidemics).

## Author

- **sechseck:** [sechseck](https://github.com/sechseck)

## Co-Author

- **fsch2:** [fsch2](https://github.com/fsch2)

## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/sechseck/epidemics/blob/master/license.md) file for details

## References

**Batista, Milan** (2020). Estimation of the final size of the coronavirus epidemic by the SIR model. [link](https://www.researchgate.net/publication/339311383_Estimation_of_the_final_size_of_the_coronavirus_epidemic_by_the_SIR_model)

https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model

## Acknowledgments

* A big Thank You goes to Milan Batista from Ljubljana for his excellent article. Thanks a lot, hvala lepa!
