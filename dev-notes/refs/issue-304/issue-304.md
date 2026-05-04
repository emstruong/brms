# brms#304 — latent variables/confirmatory factor analysis

Opened by jebyrnes on 2017-12-20T19:23:03Z — state: OPEN

---

## Original post

In a structural equation  modeling framework, we often define latent unmeasured variables with indicators.  For example

```
library(lavaan)
data(PoliticalDemocracy)

mod <- 'ind60 =~ x1 + x2 + x3'

fit <- sem(mod, data=PoliticalDemocracy)

summary(fit)
```

produces loadings and variances

```
Latent Variables:
                   Estimate  Std.Err  z-value  P(>|z|)
  ind60 =~                                            
    x1                1.000                           
    x2                2.193    0.142   15.403    0.000
    x3                1.824    0.153   11.883    0.000

Variances:
                   Estimate  Std.Err  z-value  P(>|z|)
   .x1                0.084    0.020    4.140    0.000
   .x2                0.108    0.074    1.455    0.146
   .x3                0.468    0.091    5.124    0.000
    ind60             0.446    0.087    5.135    0.000
```

Note, the loading on ind60 -> x1 is fixed to 1 for identifiability by common convention. In this model, the model is fit to the covariance matrix of the data. Think of ind60 as causing, i.e., being a predictor of, all of its indicators.

Moreover, `lavPredict` can generate predicted factor scores, for example. In full SEMs, latent variables are later used just like any other variable - as predictors, as responses, etc.

Would it be possible to incorporate factor analysis or latent variables in general into brms?

---

## Comment by paul-buerkner on 2017-12-20T20:07:09Z

I think that this should be possible. In fact I have thought about this feature for some time already. Given the complexity of such a feature (at least when combined with all the modeling flexibility of brms), it will likely take some time until I have fully figured out how a good implementation looks like. But when its done, it will put brms largely on par with lavaan and mplus when it comes to SEMs. This would really be awesome, I think.

---

## Comment by jebyrnes on 2017-12-20T20:26:12Z

Awesome! There’s a flip side as well - composites which are an estimated weighted sum of variables, but I can talk more about that once the latent bit is figured out. So great!

On Dec 20, 2017, at 3:07 PM, Paul-Christian Bürkner <notifications@github.com<mailto:notifications@github.com>> wrote:


I think that this should be possible. In fact I have thought about this feature for some time already. Given the complexity of such a feature (at least when combined with all the modeling flexibility of brms), it will likely take some time until I have fully figured out how a good implementation looks like. But when its done, it will put brms largely on par with lavaan and mplus when it comes to SEMs. This would really be awesome, I think.

—
You are receiving this because you authored the thread.
Reply to this email directly, view it on GitHub<https://github.com/paul-buerkner/brms/issues/304#issuecomment-353169104>, or mute the thread<https://github.com/notifications/unsubscribe-auth/AAuOosRZPONxih5HxGFNvJQc07Nk3Pl3ks5tCWjtgaJpZM4RI2TM>.



---

## Comment by paul-buerkner on 2017-12-20T21:04:00Z

Feel free to talk about this now. Maybe it makes sense to implement composites along with CFAs at the same time.

---

## Comment by jebyrnes on 2017-12-21T19:29:22Z

Sure - so, a composite is conceptually both simpler and odder.

Let's say you have two predictors and one response.

y ~ x1 + x2

Now, let's say instead you want a composite to represent the effect of x1 and x2 on y.

y ~ c1
c1 ~ x1 + x2

To keep the model identified, C1 has 0 residual variance. And the coefficient for y ~ c1 is fixed to 1. So, the coefficients for x1 and c2 predicting c1 really end up the same as those in the above y ~ x1 + x2 relationship.

Why would you do such a think? Well, for generality of model structure. Also, in the SEM world, then estimating the standardized coefficient between c1 and y (e.g., sd(y)/sd(c1)) is a *very* useful metric, particularly if, say, x2 = x1^2. Or if the composite influences multiple response variables, in which case the coefficients can be quite different from the simple case I presented above.

There is also a "latent composite" where one does not fix the error variance on c1 to 0, and c1 can influence multiple different endogenous variables. This is often a less table model form, but it is occasionally useful.

---

## Comment by paul-buerkner on 2017-12-21T23:17:17Z

Thanks for the explanation! I think implementing this will require very similar formula parsing than CFAs, so it will likely make sense implementing CFAs and composites together.

---

## Comment by jebyrnes on 2017-12-22T02:21:01Z

Awesome! I think the trick will be making sure that one can define a latent variable in one bf, and then use it as a predictor in another. Good luck!

---

## Comment by jebyrnes on 2019-02-19T19:39:13Z

I've been playing with this a bit more, based on the missing data vignette. I think I'm getting somewhere, and would love some opinions. The convergence is poor, and I'm wondering why and if we could improve it if it would lead to a better match. Although, notably, the variance in the factor is about the same...  But, thoughts?

```

library(brms)
library(lavaan)

model <- ' 
     ind60 =~ x1 + x2 + x3
'

fit <- sem(model, data = PoliticalDemocracy)

PoliticalDemocracy$ind60 <- as.numeric(NA)

fact_mod <- bf(x1  ~ 1*mi(ind60)) +
  bf(x2  ~ mi(ind60)) +
  bf(x3  ~ mi(ind60)) +
  bf(ind60 | mi() ~ 1)

fact_fit <- brm(fact_mod+
                  set_rescor(rescor = FALSE), data = PoliticalDemocracy, cores = 1)

summary(fit)
summary(fact_fit)
```

---

## Comment by paul-buerkner on 2019-02-20T09:19:18Z

The `1*variable` syntax does not work in brms that way (it just works in lavaan). You have to soft-fix the coefficient with a very strong prior. There is currently no way to hard fix coefficients, but this will come as part of the "actual" SEM implementation in brms 3.0. 

The following seems to work ok, but still not great. Note that the intercept of `ind60` is not identified at the same time as all the other intercepts.

```R
fact_mod <- bf(x1 ~ mi(ind60)) +
  bf(x2  ~ mi(ind60)) +
  bf(x3  ~ mi(ind60)) +
  bf(ind60 | mi() ~ 0) +
  set_rescor(rescor = FALSE)

fact_fit <- brm(
  fact_mod, data = PoliticalDemocracy,
  prior = prior(normal(1, 0.001), resp = "x1", coef = "miind60") +
    prior(normal(1, 3), resp = "x2", coef = "miind60") + 
    prior(normal(1, 3), resp = "x3", coef = "miind60"),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 2, cores = 2
)

summary(fact_fit)
```

---

## Comment by jackobailey on 2019-03-22T13:41:58Z

Could you do this using non-linear syntax instead? Then you wouldn't need to soft-fix it. Instead, you could only provide factor loadings for x2 and x3 on ind60, effectively fixing the factor loading of x1 on ind60 to 1.

---

## Comment by paul-buerkner on 2019-03-22T14:04:43Z

You are absolutely right. In any case, specifying these kinds of models in the current brms syntax is cumbersome and we definitely need a special, user-friendly syntax for it.

---

## Comment by jackobailey on 2019-03-22T21:24:27Z

I agree. My two cents, we need a generalisable syntax so that you can extend the type of latent variable modelling to categorical latent classes in the future without having to redo everything.

Given that we already specify the nature of parameters a-priori with prior distributions, it would perhaps be cool to fix parameters in the same way too.

Perhaps something like this for the political democracy example:

    fit <- bf(ind60 =~ x1 + x2 + x3)
    priors <- c(prior(normal(0, 1), class = "cfa"),
                prior(fixed(1), coef = "x1"),
                prior(normal(1, 1), class = "l"))

Where `=~` indicate a latent variable model, `fixed()` that the parameter is fixed a-priori to some value, and `l` (lambda) is the equivalent of `b` but for factor loadings, and `cfa` applies to all factors.

The nice thing about having a `fixed()` function is that it allows one to fix any element we like (e.g. fixing the variance instead of one of the factor loadings) and doesn't require that we start adding non-linear-like syntax in specifying these models (e.g. ind60 =~ 1*x1 + x2 + x3).

---

## Comment by paul-buerkner on 2019-03-23T09:03:36Z

I generally like your ideas! Thanks for this! A few comments:

* `=~` is no valid R operator but rather translates to `y = <one-sided formula>` which is reasonable, but nothing specific to latent variable models. Accordingly, we cannot use the `bf` function directly to specify a CFA but have to use a new function such as `lv` (as you suggested before editing your post).
* I am not sure about the `fixed` syntax via the prior, since, as things currently stand, `prior` will only apply to parameters which are actually estimated and I am not sure how well this generalizes to fixing values. In particular, users will expect this to work also for other (non-CFA related) model parameters and then be confused that it won't. Rather, we could do something like. `lv(ind60 ~ x1 + x2 + x3, loadings = c(1, NA, NA))` to fix factor loadings.

---

## Comment by jackobailey on 2019-03-23T13:40:44Z

Glad that you liked it!

In response, I think that fixing needs to be more flexible than just loadings. Sometimes we need to fix correlations and variances too. These can be complicated. As such, I think that fixing needs to be simpler than having users specify vectors of loadings/variance-covariance matrices if they want to fix parameters.

Most SEM users are probably already aware of either lavaan or MPlus syntax. Adopting either might provide a solution. I prefer lavaan's syntax over Mplus' as it is most similar to the existing brms syntax. For example, we might want to specify a simple bifactor model. This requires that the two factors be orthogonal.  Model specification might look as such:

    lv(fac1 ~ 1*x1 + x2 + x3 + y1 + y2 + y3,
       fac2 ~ 1*y1 + y2 + y3,
       fac1 ~~ 0*fac2)

Or, if you wanted to keep each latent variable in its own function and specify correlations using a different function:

    lv(fac1 ~ 1*x1 + x2 + x3 + y1 + y2 + y3) +
    lv(fac2 ~ 1*y1 + y2 + y3) +
    lv_corr(fac1, fac2, value = 0)

---

## Comment by paul-buerkner on 2019-03-24T10:19:54Z

The general brms approach is to use R formulas for the model specification and differentiate their meanings by the functions to which they are supplied. That is, your second proposal goes into the direction I have in mind. I am not so much concerned with mirroring the syntax of any other package but to create a syntax that is consistent and intuitive within the brms framework.

Operators like `~~` can be valid R syntax but don't imply an actual operator. Instead, `~~` implies a formula defined inside a formula, which does not make a lot of sense.

The lavaan like specification `lv(fac2 ~ 1*y1 + y2 + y3)` is not too bad, but interferes with how interactions are specified in basic R formula syntax. Perhaps, we could do something like `lv(fac2 ~ fix(y1, 1)+ y2 + y3)` and `lv(fac2 ~ name(y1, "Ly1")+ y2 + y3)` to fix and name coefficients, respectively.

---

## Comment by jackobailey on 2019-03-24T14:13:27Z

I've given this a lot of thought, and I'd like to see this become a reality, so I'm going to keep making suggestions. Apologies in advance for what is likely to be an essay!

There are many possibilities open to Bayesian CFA that are difficult/impossible with Frequentist CFA (e.g. latent interactions). Taking this and the intention to make the syntax consistent with brms syntax writ-large, perhaps we're looking at this all backwards.

It would make more sense to treat latent variables like `mi()` or `me()`. That is, to switch the ordering so that the factor is on the _right_ and not the left of the equation. This would reflect the fact that what we're dealing with is just regression where the predictor is unmeasured. Model specification might be a little longer, but it would also make things more consistent, and (in that spirit) more flexible.

For example, imagine we have 7 items (x1-7) that we think are caused by two factors plus their interaction. Imagine also that we want the factors to be correlated, and that we think there is some additional correlation between x1 and x2. We could write this as:

    bf(cbind(x1, x2, x3, x4, x5, x6, x7) ~ lv(f1)*lv(f2)) +
    corr(~ f1 + f2 + f1:f2) +
    corr(~ x1 + x2) +
    fix(path(from = f1, to = x1, value = 1),
        path(from = f2, to = x1, value = 1),
        path(from = f1:f2, to = x1, value = 1))

This has a number of benefits:

1. Specifying latent interactions or latent-observed interactions is easy here.
2. If you want to include observed variables as covariates, you just add them to the right-hand side of the equation as you would with any other regression.
3. You could add arguments to `lv()` just as you would any other function that would allow users to have the latent variable follow a desired distribution. E.g. `lv(f1, family = categorical())` for a latent class model.

Further, a generic `fix()` function, might allow you to fix other aspects of the model too. E.g.:

    fix(path(from = f1, to = x2, value = "a"),
        path(from = f2, to = x2, value = "a"),
        corr(~ x1 + x2, value = 1),
        var(f1, value = "b"),
        var(f2, value = "b"),
        var(f1:f2, value = 1))

---

## Comment by paul-buerkner on 2019-03-24T15:17:02Z

Thanks again! Much of this is in accordance with the structure I have in mind, but just realized I haven't written it down somewhere. Below is the my suggestion on how to access a latent variable once defined by some measurement formula (whose syntax is to be discussed).

* Suppose we have a latent variable `x`, then we can access `x` via `lv(x)` (or perhaps even `l(x)`) to use it as predictor. For instance `bf(y ~ lv(x)*z)` if `y` and `z` as known/observed variables. This behavior is the same as for `mi` and `me` terms.
* If `x` should be used a response variable, then `bf(x | lv() ~ ...)` will do. We can then also specify the family of the latent variable if desired. This behavior is the same as for `mi` terms.
* If `x` is not modeled as response somehwere, it will have a default model, which will likely just be `x ~ normal(0, SD_x)` with `SD_x` either fixed or estimated depending on the measurement model.

If we decided on using the syntax laid out above, the only thing that remains to be discussed is the measurement model. I agree what having the observed variables on the left and the latent variable on the right hand side is more consistent and an perhaps more flexible. In that setting, we can use all the syntax we already have in place. In order to fix a parameter value for identifiability, we could use an argument `fix` (or something like that) in `lv`.  Passing a name to `fix` would allow to set values equal if their names match.

So, for two latent variables `x` and `y` each measured by 3 variables and `x` predicting `y`, we could have

```
bf(x1 ~ lv(x, fix = 1) + bf(x2 ~ lv(x)) + bf(x3 ~ lv(x)) +
  bf(y1 ~ lv(y, fix = 1) + bf(y2 ~ lv(y)) + bf(y3 ~ lv(y)) +
  bf(y | lv() ~ lv(x))
```



---

## Comment by paul-buerkner on 2019-03-24T15:18:48Z

The only thing I don't understand in your proposal is the latent interaction part. Could you perhaps write down in a set of equations how such an interaction model you specify via

```
bf(cbind(x1, x2, x3, x4, x5, x6, x7) ~ lv(f1)*lv(f2)) +
corr(~ f1 + f2 + f1:f2) +
corr(~ x1 + x2) +
fix(path(from = f1, to = x1, value = 1),
    path(from = f2, to = x1, value = 1),
    path(from = f1:f2, to = x1, value = 1))
```

would look like mathematically?



---

## Comment by jackobailey on 2019-03-24T16:22:47Z

This looks really promising!

As for the interaction, it was just an example of flexibility that Bayesian CFA offers that's hard in Frequentist CFA. The example was overly-complicated just to be able to include some of the many requirements SEM might entail.

To be honest, I find working out the maths pretty hard and it generally takes me a long time. Instead, I've made a path diagram that might explain what I mean. As all the factors load on all of the items, it's basically an exploratory factor analysis model.

<img width="882" alt="Screenshot 2019-03-24 at 16 09 25" src="https://user-images.githubusercontent.com/37049065/54882168-34feaf80-4e4f-11e9-9683-3ee818922904.png">

_(edit: I've just realised that I modified an existing Frequentist path diagram I had for this and have left the residual errors in. Obviously, this doesn't make sense in Bayesian SEM. I'm going to leave them in though because I think my intentions are clear and it's too much effort to amend the diagram)._

Here _eta_1 x eta_2_ is the interaction between _eta_1_ and _eta_2_. In a frequentist framework, such a latent interaction would require that we estimate it based on the interactions of every combination of x1-7. Other than a potentially massive number of combinations, this has other associated difficulties. For example, where items are ordinal, multiplying observed variables makes little sense. In Bayesian CFA, where we treat the latent variable as though it is missing data, I believe that latent interactions are as simple as including the interaction on the right-hand side along the other latent variables.

Stephen Martin explains it a little bit [here](http://srmart.in/bayesian-cfa-primer/).

---

## Comment by paul-buerkner on 2019-03-24T19:51:27Z

I don't think the measurement model above makes a lot of sense. The interaction between two latent variables should not have a measurement model itself but arises simply by multiplying the two latent variables once they are measured. Or am I mistaken?

In any case, I agree that with regard to latent interactions, Bayesian SEMs have a clear advantage over frequentist SEMs.

---

## Comment by jackobailey on 2019-03-24T21:53:26Z

I could well be wrong. And it's also worth keeping this thread limited to brms SEM syntax. Nevertheless, I'll explain my reasoning then we can move it. That way, you can see whether it makes sense. And, if it does, it might be an interesting test case.

Let's assume that everyone possesses [latent authoritarian-liberal and economic left-right attitudes](https://en.wikipedia.org/wiki/Political_compass). We might wish to estimate these attitudes using observed responses to political items in electoral surveys. As such, we might fit the following model (plus priors, correlations, etc.):

![gif-1](https://user-images.githubusercontent.com/37049065/54885834-82434700-4e78-11e9-8055-4e8a45e7f633.gif)

Where AL represents latent authoritarian-liberal attitudes, LR represents latent economic left-right attitudes, y_{ij} is the observation of item _j_ for subject _i_, lambdas are factor loadings, and sigmas are item-specific standard deviations.

This model seems reasonable, but we might also expect it to fall down where the item's topic is salient to _both_ authoritarian-liberal and economic left-right attitudes. For example, climate change is an authoritarian-liberal issue because authoritarians are more likely to deny climate science and liberals more likely to embrace it. But it is also a left-right issue because of the costs and market disruption involved in transitioning to a green economy. In this case, very right-wing and very authoritarian subjects might show especially negative responses and very left-wing and very liberal subjects especially positive responses above-and-beyond the latent variables' main effects alone. As such, we might instead want to specify a model that makes this interaction explicit:

![gif](https://user-images.githubusercontent.com/37049065/54885990-95571680-4e7a-11e9-9e2e-27e1d529d2b6.gif)

In these cases, I would expect taking the interaction once the latent variables are measured is not enough. This is because doing so does not account for the fact that the two interacted to induce particular responses to the measurement items in the first place. Again, this might not make any sense and I might be mistaken, but it doesn't seem too unreasonable an expectation. I should add as a final note that not all items might load onto all factors in this case (these models would probably not be identified), but that I have simply assumed they do for notational simplicity.

---

## Comment by jebyrnes on 2019-03-25T00:35:29Z

`I would expect taking the interaction once the latent variables are measured is not enough. This is because doing so does not account for the fact that the two interacted to induce particular responses to the measurement items in the first place. `

Huh. I think I disagree here. Thinking about the definition of a latent variable, it is an unmeasured quantity that gives rise to its indicators. An interaction between two different LVs should provide no problem. One of the reasons we get into thinking about how interactions load on to their components indicators in covariance based SEM is due to all of the bookkeeping we have to do in keeping track of correlations that arise in our model predicted covariance matrix that must match those in our observed covariance matrix - and if you're using a latent interaction, it can wreak havoc. But, piecewise approaches - both bayesian and non (which we're still thinking about - although this thread has helped convince me the item that gives me pause might not be a problem - more on that in a separate post) don't have to do this bookkeeping, as they are not comparing observed and model produced covariance matrices, but rather we can look at how well our posterior distributions match those of indicators.

So, I don't think it's a problem.

---

## Comment by jebyrnes on 2019-03-25T00:39:27Z

The problem I've been worried about is latent response variables.  Consider 

Y =~ y1 + y2 + y3
Y ~ x
(to use lavaan syntax)

If you fit this using covariance based methods and then fit this using piecewise methods (estimate the latent variable factor scores, run the regression) you 

1. Get a different path coefficient for Y ~ x
2. Get different factor loadings.

This is because in covariance based estimation, you're accounting for the correlation between x and y1, y2, and y3 when you estimate all of those quantities, not just estimating each piece separately and not conditioning on that correlation.

I admit, when I first saw this (I was working on some code to do LV for `piecewiseSEM`), I put my keyboard down and walked away. And yet, this conversation is making me think that this might be a feature, not a bug, of moving from covariance to equation-level estimation.  Thoughts?

---

## Comment by jebyrnes on 2019-03-25T00:44:36Z

BTW, @paul-buerkner, to go back to the very beginning, a syntax like

`bf(cbind(x1 + x2 + x3) ~ lv(ind60, fix = x1))`

is elegant and lovely, I think.

Do you think we need to use `lv()` when `ind60` is used in other places?

Say, if it's a response variable, is `lv()` still needed? Or if it's used as a predictor in another relationship? Again, I don't totally know how things work under the hood, but if it can be avoided, it would probably save some typos and such in the future.

Last thought - as a test example (and for documentation), I def. recommend using the example in the lavaan `sem` help file, as it is a classic.

```
## The industrialization and Political Democracy Example 
## Bollen (1989), page 332
model <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'
```

---

## Comment by paul-buerkner on 2019-03-25T08:59:56Z

Let's put the interaction discussion aside for a second. We can come back to this later on once there is a working implementation of SEMs in brms. 

I agree with @jebyrnes assessment on why piecewise and covariance methods can lead to different results. I want to note though that brms will not implement a piecewise approach as all quantities are estimated within the same model. As such, I expect brms SEMs to yield similar results than lavaan, although we keep the advantages of a regression (instead of a covariance) approach in the sense that we don't need to do any book keeping of the implied covariance matrix.

I think we need to change and extend

```
bf(cbind(x1 + x2 + x3) ~ lv(ind60, fix = x1))
```

slightly. First, we will use `mvbind` instead of `cbind` (as of 2.8.0 actually). Second, `fix = x1` may not be specific enough as we need to tell the value at which we will fix the loading (it may not be 1 in all cases).
As such, a `fix` function would help. Something like

```
bf(mvbind(x1 + x2 + x3) ~ lv(ind60, fix = fix(x1 = 1))
```

This would then be parsed to something like

```
bf(x1 ~ lv(ind60, fix = 1) +
  bf(x2 ~ lv(ind60)) +
  bf(x2 ~ lv(ind60)) +
```

Further, I expect that `lv()` will be required in all formulas so that brms knows it should be a latent variable. It is possible that we can add some syntactical sugar on top of the actual syntax that removes the requirement for `lv()`, but in the syntax that is used to create the stan code `lv()` will still be required.

---

## Comment by jackobailey on 2019-03-25T09:40:32Z

I'm so excited by this. It looks great.

It would be neat to use `fix` both where fixing a loading to a value and to a label, in my opinion. Something like

    bf(mvbind(x1 + x2 + x3) ~ lv(ind60, fix = fix(x1 = 1, x2 = "a", x3 = "a"))

Rather than

    bf(mvbind(x1 + x2 + x3) ~ lv(ind60, fix = fix(x1 = 1), name = name(x2 = "a", x3 = "a"))

@jebyrnes' idea about using the industrialization and political democracy example is a good one.

In brms this might look like:

    bf(mvbind(x1 + x2 + x3) ~ lv(ind60, fix = fix(x1 = 1))) +
    bf(mvbind(y1 + y2 + y3 + y4) ~ lv(dem60, fix = fix(y1 = 1, y2 = "a", y3 = "b", y4 = "c"))) +
    bf(mvbind(y5 + y6 + y7 + y8) ~ lv(dem65, fix = fix(y5 = 1, y6 = "a", y7 = "b", y8 = "c"))) +
    bf(dem60 | lv() ~ lv(ind60)) + 
    bf(dem65 | lv() ~ lv(ind60) + lv(dem60))

Given how much code this is, it's probably also worth considering how to cut it down. @paul-buerkner's original idea to have `fix` be based on the order the variables appear in `mvbind` would work. That way, `lv(dem60, fix = fix(y1 = 1, y2 = "a", y3 = "b", y4 = "c")),` which is long, could be written as `lv(dem60, fix = c(1, "a", "b", "c"))`. Further, in an extreme case where we had a long list of _j_ items and wanted to fix only the first and last loadings, we could write `fix = c(1, rep(NA, j-2), "a")`. Updating the syntax above gives:

    bf(mvbind(x1 + x2 + x3) ~ lv(ind60, fix = c(1, NA, NA))) +
    bf(mvbind(y1 + y2 + y3 + y4) ~ lv(dem60, fix = c(1, "a", "b", "c"))) +
    bf(mvbind(y5 + y6 + y7 + y8) ~ lv(dem65, fix = c(1, "a", "b", "c"))) +
    bf(dem60 | lv() ~ lv(ind60)) + 
    bf(dem65 | lv() ~ lv(ind60) + lv(dem60))

Which looks clearer. Another approach (though I don't know how simple this would be to implement) would be to borrow from MPlus and use the @-symbol to specify a fix and () to specify a label (or any other arbitrary/available symbols). This is appealing as users are already writing out all of the observed variables anyway and would not have to deal with providing NAs. Something like:

    bf(mvbind(x1@1+ x2 + x3) ~ lv(ind60)) +
    bf(mvbind(y1@1 + y2("a") + y3("b") + y4("c")) ~ lv(dem60)) +
    bf(mvbind(y5@1 + y6("a") + y7("b") + y8("c")) ~ lv(dem65)) +
    bf(dem60 | lv() ~ lv(ind60)) + 
    bf(dem65 | lv() ~ lv(ind60) + lv(dem60))

Though this might present problems where both latent and observed variables or multiple latent variables affect a given item. 

Whatever decision you make, it leaves only working out how to specify correlations between observed variables and latent factors so that we can include the following:

    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8


---

## Comment by paul-buerkner on 2019-03-26T07:50:55Z

I think users should be able to identify loadings both by name and by position, that is both `fix(x1 = 1, x3 = "a")` and `fix(1, NA, "a")` should work. I don't mind if the syntax is a little bit more verbose than that of lavaan or Mplus as long as it is easy to remember and consistent. 

The Mplus syntax via `@1` and `("a")` is not only R-like, but also has a problem (as you notice) when multiple latent variables loading on the same item. In other words, in the regression formulation (manifest left, latent right), we are currently advocating, the loading has to be a property of the latent variable. If it is the other way round, the loading is a property of the manifest variable instead.

Let's discuss the next set of questions:

* How do specify (residual) correlations between manifest and/or latent variables?
* How to specify grouping, i.e. when we expect different groups to have different factor loadings?
* How to specify clustering for multilevel SEMs?
* How to do predictions for new data that involve computing estimates of the latent variable first?

Suggestions are welcome!

---

## Comment by jackobailey on 2019-03-26T08:36:26Z

I have suggestions for 1 and 2.

1. How about keeping it simple and having a general function to specify correlations between variables: e.g. `cor(x, y)`, `cor(lv(x), lv(y))`. Where we want to specify correlations between many variables (as might be the case with longitudinal data), we could write `cor(y1, y2, y3, y4)`. Alternatively, you might prefer a formula-like syntax `~ y1 + y2 + y3 + y4` for the same purpose. We could even include a `fix` argument so that we can fix/label correlations. E.g. `cor(x, y, fix = 1)`.

2. For grouping, what about `bf(mvbind(x1, x2, x3, x4) ~ lv(lat, group = group, fix = fix(x1 = 1, x2 = c("a", "b"), x3 = c("a", NA))))`? That way we could specify a lot of these group-specific loadings as objects before placing them in `lv`.

---

## Comment by bwiernik on 2019-05-10T14:50:30Z

The lavaan syntax is pretty well-embedded in R SEM circles and is also being supported by the other major R SEM package, OpenMx, so it would be really nice if brms at least supported the option of providing an SEM model using lavaan syntax. 

Cf https://github.com/tbates/umx/issues/24#issuecomment-491312130

---

## Comment by paul-buerkner on 2019-05-10T15:59:39Z

It could indeed be worth it to add lavaan syntax (or parts of it) as syntactical sugar on top. However,
my plan is to support a lot of SEM models which are out of scope of lavaan (or Mplus) and as such it cannot be the primary syntax for SEMs in brms.

@jackobailey Sorry for not responding to your ideas. I will have to think of the syntax for correlations more thoroughly, but you ideas help me a lot in this regard!

---

## Comment by bwiernik on 2019-05-10T16:43:10Z

An approach would be to have a function that takes a lavaan string and returns the corresponding brms syntax. 

---

## Comment by paul-buerkner on 2019-05-10T17:00:45Z

yes. that's what I have in mind as well.

Brenton M. Wiernik <notifications@github.com> schrieb am Fr., 10. Mai 2019,
18:43:

> An approach would be to have a function that takes a lavaan string and
> returns the corresponding brms syntax.
>
> —
> You are receiving this because you were mentioned.
> Reply to this email directly, view it on GitHub
> <https://github.com/paul-buerkner/brms/issues/304#issuecomment-491354073>,
> or mute the thread
> <https://github.com/notifications/unsubscribe-auth/ADCW2AAWI6XQH2C2VBJ43HDPUWQZ7ANCNFSM4EJDMTGA>
> .
>


---

## Comment by JCruk on 2021-01-14T17:56:28Z

Is this effort still in scope?
I'm working on a factor analysis for ordinal data and have been guided by this discussion and some [other posts](https://scottclaessens.github.io/blog/2020/brmsLV/) to attempt _hacking_ brms for the task at hand.

Is there an alpha branch or anything I could test to assist in SEM development efforts?

---

## Comment by paul-buerkner on 2021-01-15T14:52:51Z

It is still in scope. Unfortunately, the grant proposal that would have involved implementing this in brms was recently rejected and I will wait until the resubmission is (hopefully) accepted.

---

## Comment by JCruk on 2021-01-15T18:00:03Z

Thanks, @paul-buerkner . Sorry to hear the grant wasn't awarded!

I'm not sure if this discussion is better suited to discourse, but I thought I'd start here to see if what I'm thinking is at all reasonable.

I'm planning some psychometric evaluation of a questionnaire where responses are on the same 7-point ordinal scale (from strongly disagree to strongly agree).

Rather than treat these data as metric, I want to try using an ordered probit model with the following _hack_ based on what I've seen:

```
# add empty latent variable to data frame
d$X <- as.numeric(NA)

KPES_Data$F1 <- cbind(KPES_Data$Q1, KPES_Data$Q2)
KPES_Data$F2 = cbind(KPES_Data$Q3,KPES_Data$Q4, KPES_Data$Q5)
KPES_Data$F3 = cbind(KPES_Data$Q6,KPES_Data$Q7, KPES_Data$Q8,KPES_Data$Q9 )

fact_mod  <- bf(F1  ~ 1 + mi(X) + (1|q|Subject)) +# factor loading 1 (constant)
        bf(F2  ~ 1 + mi(X) + (1|q|Subject)) + # factor loading 2 
        bf(F3  ~ 1 + mi(X) + (1|q|Subject)) + # factor loading 2
        bf(X | mi() ~ 1) +
        set_rescor(rescor = FALSE)

fit <- brm(fact_mod, data = KPES_Data, 
           prior = c(prior(normal(1, 0.000001), coef = miX, resp = F1),
                     prior(normal(1, 1), coef = miX, resp = F2),
                      prior(normal(1, 1), coef = miX, resp = F2)),
           family = "cumulative"(probit, threshold = "flexible"),
           iter = 2000, warmup = 1000, chains = 2, cores = 2,
           control = list(adapt_delta = 0.99, max_treedepth = 12))

```

Does that seem like a somewhat sound approach to a latent factor analysis with the current version of brms? I've set the priors on the response (question) that is thought to best represent the latent factor.


---

## Comment by JCruk on 2021-01-16T02:42:38Z

Unfortunately that didn't work because `mi` is not supported for cumulative models.

Any suggestions on how a latent model could work with ordinal data?

---

## Comment by paul-buerkner on 2021-01-16T09:34:50Z

Please ask this question on discourse. For now it probably won't work with brms.

---

## Comment by bwiernik on 2021-03-14T16:45:50Z

Just wanted to note that psychologists will want to be able to compute approximate fit indices such as RMSEA, CFI, SRMR, etc. if brms is to be adopted for SEM-type analyses. Some recent theoretical work on those indices in a Bayesian framework are https://doi.org/10.1177/0013164417709314 and https://doi.org/10.1080/10705511.2020.1764360

---

## Comment by jackobailey on 2021-06-24T09:33:41Z

Just a heads up: I think I worked out a way to fit rudimentary confirmatory factor analysis models in brms using the constant() and mi() functions. Keen on any peer review too (I might have messed up but seems to work for me). https://discourse.mc-stan.org/t/confirmatory-factor-analysis-using-brms/23139

I've also extended it to work with ordered outcomes too and can share the code if it's helpful.

---

## Comment by paul-buerkner on 2021-06-25T15:39:38Z

The brms code you provide looks reasonable to me. I understand it provides the expected results?

---

## Comment by jackobailey on 2021-06-25T16:05:45Z

Yeah, it does. Output below. Picked bad effect sizes, so they're scaled by .5 here due to the constraints on the factor loadings.

```
 Family: MV(gaussian, gaussian, gaussian, gaussian) 
  Links: mu = identity; sigma = identity
         mu = identity; sigma = identity
         mu = identity; sigma = identity
         mu = identity; sigma = identity 
Formula: y1 ~ 0 + mi(xo) 
         y2 ~ 0 + mi(xo) 
         y3 ~ 0 + mi(xo) 
         xo | mi() ~ 1 
   Data: dta (Number of observations: 200) 
Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
         total post-warmup samples = 4000

Population-Level Effects: 
             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
xo_Intercept     0.16      0.15    -0.14     0.46 1.00     6213     3108
y2_mixo          0.51      0.05     0.42     0.60 1.00     4023     3189
y3_mixo          0.26      0.04     0.18     0.34 1.00     6192     3420
y1_mixo          1.00      0.00     1.00     1.00 1.00     4000     4000

Family Specific Parameters: 
         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sigma_xo     1.97      0.12     1.74     2.21 1.00     4100     3358
sigma_y1     1.00      0.00     1.00     1.00 1.00     4000     4000
sigma_y2     1.00      0.00     1.00     1.00 1.00     4000     4000
sigma_y3     1.00      0.00     1.00     1.00 1.00     4000     4000

Samples were drawn using sample(hmc). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

---

## Comment by jebyrnes on 2021-12-06T16:43:08Z

Hey, all, I fell out of this for a while due to other issues in professional life, but I wanted to say how amazing this discussion is.

1. https://discourse.mc-stan.org/t/confirmatory-factor-analysis-using-brms/23139 is amazing @jackobailey - I'm going to use it to teach with, and I see no reason why it shouldn't be a foundation for future work!

2. The `cor(y1, y2, y3, y4)` syntax seems ideal - and could address something I've been thinking about with respect to `set_rescor()` - which dovetails with https://github.com/paul-buerkner/brms/issues/957 - the need for more flexible covariance matrices to accomodate things like correlated error between some terms, but not others. Which is a real need for some models!

3. @paul-buerkner - bummer about not getting this funded! It should be! If there is anything I can do to help in terms of providing a biological motivation, a teaching motivation (since I teach this a lot!), or any source of collaboration that would help you get funded, I am more than happy to help! We're trying to grapple with this in frequentist piecewiseSEM, and Shipley and Douma's recent paper on MAGs as a way of shoehorning latent variables and correlated errors into this framework was really really great - https://www.tandfonline.com/doi/abs/10.1080/10705511.2020.1871355 - but, focused on testing and honestly makes the analyst have to do some mental backflips that seem unnecessary - and are in a Bayesian Framework where we can build the model we mean instead of trying to figure out how to translate between DAGs and MAGs.

I guess all of which to say is, I'm cheerleading and let me know if I can help at all here - both in terms of providing code or providing anything in the none-code world we live in. HA!

---

## Comment by andrjohns on 2022-03-17T04:46:37Z

I've also been thinking that once I've finished up the Gaussian copulas from #1317, that it would make a natural base for implementing LV models (in the generated Stan code, the `brms` formula side is a bigger can of worms!). Given that the copula approach allows for pretty flexible modelling of correlated outcomes, we could implement latent variable models that work with any combination of input variable types (as long as they have a CDF).

So if we use the marginal likelihood parameterisation for the model:
![image](https://user-images.githubusercontent.com/27717896/158738084-7b4e86d2-9206-49ae-98fc-d444f6c38b9e.png)

Then we can specify the model using a multivariate copula with outcome correlations defined by Sigma and the intercepts/variances (where relevant) handled by the respective marginal distributions. This means that `brms` could have a single framework for latent variable models that can then be applied to all outcome types

Bayesian LV's with ordinal outcomes would be a killer feature for `brms` for some of my own work, so I think this provides a pretty exciting avenue for getting that implemented!

---

## Comment by paul-buerkner on 2022-03-20T18:22:05Z

I agree. I have recently aquired a grant which involves latent variable models within brms so copulas will be perfect for that. In fact I actually mention copulas in the proposal so your timing on this is just perfect :-D

Edit: Just to clarify, I do not intent to use a covariance formulation for latent variable models but rather a conditional formulation via explicit latent variables, because I perceive the latter as more flexible and easier to integrate with all the other brms features. Still, for all the residual correlations, we would need copulas.

---

## Comment by andrjohns on 2022-03-21T07:23:01Z

Oh excellent news, thats great to hear! Yeah the conditional formulation should work just as easily with the latent variables being used by the marginal distributions. 

Out of interest, would the conditional parameterisation work with cross-validation? Given that the individual latent variable score can't be passed as an input for the prediction, is it just as equivalent to first sample from a normal RNG and then use that for the individual score, or would the marginal parameterisation be used for the cross-validation side of things instead?

---

## Comment by paul-buerkner on 2022-03-21T07:45:19Z

That is a good question. We could hope that PSIS does the work on finding the correct LOO-posterior even if in that LOO-posterior the individual latent variable score didn't have any data to back it up. Perhaps using a marginal parameterization could make sense there but I really don't know yet for sure.

The reason I don't like the marginal parameterization is that is, to my current understanding, not super easy to combine with all the other features of brms that work on the mean rather than the covariance (e.g. splines, GP, multilevel structure etc.) and also doesn't play well necessarily with distributional regression, which is part of the project as well. So I think the conditional formulation is the way to go despite creating some potential problems for CV.

---

## Comment by ecmerkle on 2022-03-22T14:18:57Z

I just stumbled on this issue and wanted to say that I would be interested in helping with the lavaan-to-brms translation and could provide lots of lavaan test models.

Also, I have looked at CV for these models and think it is almost more difficult than the model estimation itself. The marginal likelihood seems preferable (see [this paper](https://arxiv.org/abs/1802.04452)) but is tough to approximate for non-normal models. I have had some success with quadrature and with post-estimation importance sampling, but they are generally slow and maybe difficult to get working with all the models that brms could support. The copula idea [above](https://github.com/paul-buerkner/brms/issues/304#issuecomment-1070322293) might be worthwhile for this, though I am not sure that all the possible brms models can be covered by a lisrel-like framework.

---

## Comment by paul-buerkner on 2022-03-23T07:32:12Z

A lavaan-to-brms translation would be really nice, thank you! I am currently brainstorming in my head a bit how the brms-native syntax shall look like and then we should be able to add syntactical sugar on top of this, for example, a translation from (some parts of) lavaan syntax to the brms models. I will let you know once I have made progress on the former. 

I agree with regard to CV. Sounds like there is still a lot of work to do to make CV usable for latent variable models.

---

## Comment by E-M-McCormick on 2022-12-01T22:10:22Z

Got sent this thread randomly, so you might be well past needing any additional help with this sort of effort, but if you ever do need help with anything related to this implementation, lavaan translation, or just need a monkey to try to test/break things, please let me know. Moving out from under Mplus for models currently unavailable in lavaan would be a godsend.

---

## Comment by paul-buerkner on 2022-12-02T12:34:14Z

thank you! I will let you know if we need additional help on this!

Ethan M McCormick ***@***.***> schrieb am Do., 1. Dez. 2022,
23:10:

> Got sent this thread randomly, so you might be well past needing any
> additional help with this sort of effort, but if you ever do need help with
> anything related to this implementation, lavaan translation, or just need a
> monkey to try to test/break things, please let me know. Moving out from
> under Mplus for models currently unavailable in lavaan would be a godsend.
>
> —
> Reply to this email directly, view it on GitHub
> <https://github.com/paul-buerkner/brms/issues/304#issuecomment-1334512585>,
> or unsubscribe
> <https://github.com/notifications/unsubscribe-auth/ADCW2AFJ7YFVCJVHVNMNDB3WLEO5VANCNFSM4EJDMTGA>
> .
> You are receiving this because you were mentioned.Message ID:
> ***@***.***>
>


---

## Comment by GidonFrischkorn on 2024-04-04T12:47:41Z

Like others before me, I just read through the whole thread and think it would be awesome to have SEM/CFA features in `brms`. One thing I was not sure if it was discussed so far is to combine parameters from a regression model and their variances/covariances into a SEM model. 

More specifically, I mean providing structural equations for the covariance (G) matrix  of estimated parameters of a GLM and thus calculate CFAs/SEMs on these. Ideally, it would be possible to include measurement models for external covariates to estimate relationships of these with the factors extracted from the model parameters. In my mind, this goes in the direction of Multi-Level-SEM, where parts of the model can be a GLM and parts of the model are a SEM. Maybe, that is what you meant @paul-buerkner above, when discussing the marginal parametrization. If that's the case then you will likely already work on this.

One challenge I anticipate in this is that not all variables in the model need to stem from the same underlying data distribution. For example, some responses in such models could be accuracies (with binomial distributions), while others are likert-scale responses (cumulative), and yet others reaction times (log-normal or inv.gaussian).

In any case, I would also be happy to contribute to this feature, if help is needed, as I would love to see this in `brms`.

---

## Comment by karchjd on 2026-02-20T11:10:59Z

It would be awesome to have SEM features in brms, freeing users from reliance on Mplus for complicated models. Has there been any progress on this? Are there plans or a rough timeline?

I know it’s a big project. Adding latent variables while keeping all the existing brms features sounds very challenging. Just curious where things stand.


---

## Comment by paul-buerkner on 2026-02-20T11:19:44Z

In some form yes. You can use `re()` terms on the brms3 branch to use random effects as predictors in other parts of the model. We will have a paper on this soon. Not super convenient for all kinds of SEMs but together with using 'mi()' provides quite a bit of flexibility even if the syntax is not yet super crisp.

---

## Comment by karchjd on 2026-02-20T11:22:43Z

Thanks for the swift response! That's great news! I assume a draft of the paper is not yet available as a preprint?

---

## Comment by paul-buerkner on 2026-02-20T11:47:44Z

Nope, but in a few weeks.

