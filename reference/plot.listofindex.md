# Plot community assembly index

Plot community assembly index and confidence intervals using a list of
index. S3 method for class listofindex.

## Usage

``` r
# S3 method for class 'listofindex'
plot(x, type = "normal", 
  col.index = c("red", "purple", "olivedrab3"), add.conf = TRUE, 
  color.cond = TRUE, val.quant = c(0.025, 0.975), 
  grid.v = TRUE, grid.h = TRUE, xlim = NULL, ylim = NULL, 
  cex.text = 0.8, plot.ask = FALSE, srt.text = 90, alpha = 0.4, ...)
```

## Arguments

- x:

  A list of index and related null models obtained from to the
  as.listofindex function.

- type:

  Type of plot. Possible type = "simple", "simple_range", "normal",
  "barplot" and "bytraits".

- col.index:

  Vector of colors for index.

- add.conf:

  Logical value; Add confidence intervals or not.

- color.cond:

  Logical value; If color.cond = TRUE, color points indicate
  T-statistics values significatively different from the null model and
  grey points are not different from null model.

- val.quant:

  Numeric vectors of length 2, giving the quantile to calculate
  confidence interval. By default val.quant = c(0.025,0.975) for a
  bilateral test with alpha = 5%.

- grid.v:

  Logical value; print vertical grid or not

- grid.h:

  Logical value; print horizontal grid or not

- xlim:

  Numeric vectors of length 2, giving the x coordinates range

- ylim:

  Numeric vectors of length 2, giving the y coordinates range

- cex.text:

  Numeric value; the magnification to be used for text relative to the
  current setting of cex

- plot.ask:

  Logical value; ask for plotting the next plot or not.

- srt.text:

  Degree of rotation for text.

- alpha:

  Degree of transparency for null models aera.

- ...:

  Any additional arguments are passed to the plot function creating the
  core of the plot and can be used to adjust the look of resulting
  graph.

## Details

S3 method plot for class listofindex: -Normal type plot means, standard
deviations, ranges and confidence intervals of T-statistics.
-Simple_range type plot means, standard deviations and range of
T-statistics -Simple type plot T-statistics for each site and traits and
the mean confidence intervals by traits -Barplot type plot means,
standard deviations and confidence intervals of T-statistics in a
barplot fashion -Bysites type plot each metrics for each sites -Bytraits
type plot each metrics for each traits

## Value

None; used for the side-effect of producing a plot.

## Author

Adrien Taudiere

## See also

[`as.listofindex`](https://adrientaudiere.github.io/cati/reference/as.listofindex.md);
[`plot.Tstats`](https://adrientaudiere.github.io/cati/reference/Tstats.md);
[`ses.listofindex`](https://adrientaudiere.github.io/cati/reference/ses.listofindex.md)

## Examples

``` r
  data(finch.ind)

  res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
  sp = sp.finch, nperm = 9, print = FALSE)
#> Warning: This function exclude 1137 Na values

  if (FALSE) { # \dontrun{
    #### Use a different regional pool than the binding of studied communities
    #create a random regional pool for the example
  
    reg.p <- rbind(traits.finch, traits.finch[sample(1:2000,300), ])
  
    res.finch2 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
      sp = sp.finch, reg.pool=reg.p, nperm = 9, print = FALSE)  
      
      plot(as.listofindex(list(res.finch,res.finch2)))
    
    
      #### Use a different regional pool for each communities
    #create a random regional pool for each communities for the example
    list.reg.p <- list(
    traits.finch[sample(1:290,200), ], traits.finch[sample(100:1200,300), ], 
    traits.finch[sample(100:1500, 1000), ], traits.finch[sample(300:800,300), ],
    traits.finch[sample(1000:2000, 500), ], traits.finch[sample(100:900, 700), ] )

    # Warning: the regional pool need to be larger than the observed communities
    table(ind.plot.finch)
    # For exemple, the third community need a regional pool of more than 981 individuals
    
    res.finch3 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
      sp = sp.finch, reg.pool=list.reg.p, nperm = 9, print = FALSE)  
      
      plot(as.listofindex(list(res.finch, res.finch2, res.finch3)))  
  } # }
```
