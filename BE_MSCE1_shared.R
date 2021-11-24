#### time increments
s = .5; tau = s/2; time.increment = s

#### pre-define p0 and alphaXsi for function genSize
xsi = (exp(-p*s) - exp(-q*s))/((q+alpha)*exp(-p*s)-(p+alpha)*exp(-q*s))
alphaXsi = alpha*xsi
top = xsi*(alpha+p)*(alpha+q)*(q*exp(-p*s)-p*exp(-q*s))
bot = q*(alpha+p)*exp(-p*s) - p*(alpha+q)*exp(-q*s)
p0 = top/bot

#### pre-define unconditional p0 when (mu=0) and alphaXsi for function genSize
p.marginal = alpha.I*(gamma.I-1)
alphaXsi.marginal = (exp(-p.marginal*s)-1)/(exp(-p.marginal*s)-gamma.I)
xsi.marginal = alphaXsi.marginal/alpha.I
p0.marginal = xsi.marginal*(alpha.I+p.marginal)

#### p0.tau: probability of extinction of single cell after time tau
p.marginal = alpha.I*(gamma.I-1)
alphaXsi.tau = (exp(-p.marginal*tau)-1)/(exp(-p.marginal*tau)-gamma.I)
xsi.tau = alphaXsi.tau/alpha.I
p0.tau = xsi.tau*(alpha.I+p.marginal)

#### for malignant cells (marginal on detection)
p.M = alpha.M*(gamma.M-1)
alphaXsi.M = (exp(-p.M*s)-1)/(exp(-p.M*s)-gamma.M)
xsi.M = alphaXsi.M/alpha.M
p0.M = xsi.M*(alpha.M+p.M)

#### Biopsy detection 
cells.per.cm2 = X/(3*7.5)  # no. of cells per square cm (3cm BE length, 7.5cm circ.)
# area.per.biopsy = 3.75  # 4 Q biopsies every 2 cm [Area per biopsy in BE segment (which is not the biopsy area)]
area.per.biopsy =   15    # 1 Q biopsy every 2 cm
scaleDetection = cells.per.cm2*area.per.biopsy 

#### if preferred, put all parameters in a data frame 
# varSet = data.frame(increment=increment,
#                     theta=theta,
#                     p=p,q=q,alpha=alpha,mu=mu,
#                     p0=p0,alphaXsi=alphaXsi,
#                     p0.marginal=p0.marginal,alphaXsi.marginal=alphaXsi.marginal,
#                     p0.M=p0.M,alphaXsi.M=alphaXsi.M)  


