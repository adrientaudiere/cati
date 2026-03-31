# Transform index results in a list of index

Transform various results from functions Tstats, ComIndex or
ComIndexMulti in a list of index. Useful to use the functions
plot.listofindex (S3 method) and ses.listofindex.

## Usage

``` r
as.listofindex(x, namesindex = NULL)
```

## Arguments

- x:

  A list of objects of class Tstats, ComIndex or ComIndexMulti

- namesindex:

  Optionnal, the names of index in the same order as in x.

## Value

A list of observed values and corresponding "null" values (i.e. produced
by null models) in the form "list(index1, null model index1, index2,
null model index2 ...)"

## Author

Adrien Taudiere

## See also

[`ses.listofindex`](https://adrientaudiere.github.io/cati/reference/ses.listofindex.md);
[`plot.listofindex`](https://adrientaudiere.github.io/cati/reference/plot.listofindex.md)

## Examples

``` r
  data(finch.ind)

  res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
  sp = sp.finch, nperm = 9, print = FALSE)
#> Warning: This function exclude 1137 Na values

  # \donttest{
    #### Use a different regional pool than the binding of studied communities
    #create a random regional pool for the example
  
    reg.p <- rbind(traits.finch, traits.finch[sample(1:2000,300), ])
  
    res.finch2 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
      sp = sp.finch, reg.pool=reg.p, nperm = 9, print = FALSE)  
#> Warning: This function exclude 1137 Na values
      
      plot(as.listofindex(list(res.finch,res.finch2)))

    
    
      #### Use a different regional pool for each communities
    #create a random regional pool for each communities for the example
    list.reg.p <- list(
    traits.finch[sample(1:290,200), ], traits.finch[sample(100:1200,300), ], 
    traits.finch[sample(100:1500, 1000), ], traits.finch[sample(300:800,300), ],
    traits.finch[sample(1000:2000, 500), ], traits.finch[sample(100:900, 700), ] )

    # Warning: the regional pool need to be larger than the observed communities
    table(ind.plot.finch)
#> ind.plot.finch
#>     DMaj    EspHd FlorChrl  GnovTwr MrchBndl SCruInde 
#>       50      267      981      258      270      687 
    # For exemple, the third community need a regional pool of more than 981 individuals
    
    res.finch3 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
      sp = sp.finch, reg.pool=list.reg.p, nperm = 9, print = FALSE)  
#> Warning: This function exclude 1137 Na values
      
      plot(as.listofindex(list(res.finch, res.finch2, res.finch3)))  

  # }

```
