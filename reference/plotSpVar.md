# Plot populations values against species values

Plot populations values against species values. The objectif is to see
the contribution of intra-specific vs inter-specific variation to trait
gradient.

## Usage

``` r
plotSpVar(traits = NULL, ind.plot = NULL, sp = NULL, variable = NULL, 
  col.ind = rgb(0.5, 0.5, 0.5, 0.5), col.pop = NULL, col.sp = NULL, 
  col.site = NULL, resume =  FALSE, p.val = 0.05, min.ind.signif = 10, 
  multipanel = TRUE, col.nonsignif.lm = rgb(0, 0, 0, 0.5), 
  col.signif.lm = rgb(1, 0.1, 0.1, 0.8), silent =  FALSE)
```

## Arguments

- traits:

  Individual Matrix of traits with traits in columns.

- ind.plot:

  Factor defining the name of the plot in which the individual is.

- sp:

  Factor defining the species which the individual belong to.

- variable:

  A matrix of variables corresponding to each site (in rows) and each
  trait (in columns). If you want to plot all traits against one
  variable, variable can be a vector of numerical values.

- col.ind:

  Color for individual values.

- col.pop:

  Color for populational mean values.

- col.sp:

  Color for species mean values.

- col.site:

  Color for sites mean values.

- resume:

  Logical, if TRUE plot a simple form of the plot.

- p.val:

  Choosen p.value to print significant linear relationship using linear
  model. Argument past to the lm funtion internally.

- min.ind.signif:

  Minimum individual to print significant linear relationship.

- multipanel:

  Logical value. If TRUE divides the device to shown several traits
  graphics in the same device.

- col.nonsignif.lm:

  Color for non significant linear relationship.

- col.signif.lm:

  Color for significant linear relationship.

- silent:

  Logical value, if resume = FALSE do not print warning argument.

## Value

None; used for the side-effect of producing a plot.

## Author

Adrien Taudiere

## See also

[`plotDistri`](https://adrientaudiere.github.io/cati/reference/plotDistri.md)

## Examples

``` r
  data(finch.ind)
  
  #Random variable for this example
  variable <- c(1,5,15,6,3,25)
  # \donttest{
  plotSpVar(traits.finch, ind.plot.finch, sp.finch, variable, 
  silent = TRUE)


  #If we change the value of the threshold 
  #(alpha = 10% instead of 5% 
  #and the minimum individual to represent significativity 
  #fixed to 3 instead of 10 by default) 
  #we can see some significant relationships.

  plotSpVar(traits.finch, ind.plot.finch, sp.finch, variable, 
  p.val = 0.1,  min.ind.signif = 3, silent = TRUE)



  #For a more simple figure, add the option resume = TRUE. 
  #Again if we change the value of the threshold 
  #(alpha = 10% instead of 5% 
  #and the minimum individual to represent significativity
  # fixed to 3 instead of 10 by default) 
  #we can see some significant relationships.

  plotSpVar(traits.finch, ind.plot.finch, sp.finch, variable, 
  silent = TRUE, resume = TRUE, col.pop = "grey")

  
  plotSpVar(traits.finch, ind.plot.finch, sp.finch, variable, 
  silent = TRUE, resume = TRUE, col.pop = "grey", col.sp = "black")

  
  plotSpVar(traits.finch, ind.plot.finch, sp.finch, variable, 
  silent = TRUE, resume = TRUE, col.pop = "grey", col.sp = "black", 
  p.val = 0.1,  min.ind.signif = 3)

# }
```
