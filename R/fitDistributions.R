# fits the double geometric distribution
fitDoubleGeometric<-function(n,trials=1000,min=NA,max=NA)	{
	n = sort(n,decreasing=T)
	spp = length(n)
	if (is.na(min))
		min = spp
	if (is.na(max))
		max = 10 * spp
	p = n / sum(n)
	best = 999999
	k = 0.5
	bestk = 1
	r = min
	bestr = min - 1
	while (r <= max && r == bestr + 1)	{
		for (z in 1:trials)	{
			lastk = k
			k = k + rnorm(1,sd=1 / z)
			if (k <= 0)
				k = lastk
			dg = array()
			dg[1] = k
			for (i in 2:r)	{
				if ( i > r - i + 1 )	{
					dg[i] = dg[i-1] * k ** (1 / sqrt( (r - i + 1) / r ))
				} else	{
					dg[i] = dg[i-1] * k ** (1 / sqrt( i / r ))
				}
			}
raw = dg
			dg = dg / sum(dg)
			kl = sum(p * log(p / dg[1:spp]))
			if (kl < best)	{
				best = kl
				bestk = k
				bestr = r
				bestdg = dg
			}
			if (bestk != k)	{
				k = lastk
			}
		}
		r = r + 1
	}
	return(list(distribution=bestdg,richness=bestr,k=bestk,fit=best))
}

# fits the log normal distribution
fitLogNormal<-function(n,trials=1000,min=NA,max=NA)	{
	n = sort(n,decreasing=T)
	spp = length(n)
	if (is.na(min))
		min = spp
	if (is.na(max))
		max = 10 * spp
	p = n / sum(n)
	best = 999999
	s = 2
	r = min
	bestr = min - 1
	while (r <= max && r == bestr + 1)	{
		for (z in 1:trials)	{
			lasts = s
			s = s + rnorm(1,sd=1 / z)
			if (s <= 0)	{
				s = lasts
			}
			np = exp(sort(qnorm(seq(0.5 / r,1 - 0.5 / r,1 / r),sd=s),decreasing=T))
			np = np / sum(np)
			kl = sum(p * log(p / np[1:spp]))
			if (kl < best)	{
				best = kl
				bests = s
				bestr = r
				bestln = np
			}
			if (bests != s)	{
				s = lasts
			}
		}
		r = r + 1
	}
	return(list(distribution=bestln,richness=bestr,sd=bests,fit=best))
}

# fits the geometric series distribution
fitGeometric<-function(n)	{
	n = sort(n,decreasing=T)
	spp = length(n)
	p = n / sum(n)
	k = exp(lm( log(n) ~ c(1:spp),weights=sqrt(n) )$coef[2])
	best = 999999
	for (i in 1:100)	{
		lastk = k
		k = k + rnorm(1,sd=0.01)
		q = array()
		for (i in 1:(10*spp))	{
			if (i > 1)
				q[i] = q[i-1] * k
			else
				q[i] = k
		}
		q = q / sum(q)
		kl = sum(p * log(p / q[1:spp]))
		if (kl < best)	{
			best = kl
			bestk = k
			bestgs = q[1:spp]
		}
		if (bestk != k)	{
			k = lastk
		}
	}
	return(list(distribution=bestgs,k=as.numeric(bestk),fit=best))
}

# fits the log series distribution
fitLogSeries<-function(n)	{
	n = sort(n,decreasing=T)
	spp = length(n)
	N = sum(n)
	p = n / N
	if (N > 10000)
		return(list(distribution=NA,alpha=NA,fit=NA))
	z = logSeriesParams(n)
	alpha = z[1]
	lsx = z[2]
	best = 999999
	besta = alpha
	for (i in 1:1000)	{
		lasta = alpha
		alpha = alpha + rnorm(1,sd=1)
		if (alpha <= 0)
			alpha = lasta
		lsx = N / (N + alpha)
		s = alpha * lsx^(1:(3 * N))/(1:(3 * N))
		cd = round(cumsum(s))
		z = 0
		q = array()
		for (j in 1:cd[1])	{
			z = z + 1
			q[z] = 1
		}
		for (breakpt in which(diff(cd) > 0))
			for (j in cd[breakpt]:(cd[breakpt+1]-1))	{
				z = z + 1
				q[z] = breakpt + 1
			}
		if (length(q) < spp)	{
			alpha = lasta
			next
		}
		q = sort(q,decreasing=T)
		q = q / sum(q)
		kl = sum(p * log(p / q[1:spp]))
		if (kl < best)	{
			best = kl
			besta = alpha
			bestls = q[1:spp]
		}
		if (besta != alpha)	{
			alpha = lasta
		}
	}
	return(list(distribution=bestls,alpha=besta,fit=best))
}

# needed by fitLogSeries
logSeriesParams<-function(n)	{
	a = 10
	olda = 0
	N = sum(n)
	z = 0
	while (abs(a - olda) > 0.0000001 && z < 1000)	{
		z = z + 1
		olda = a
		a = length(n) / log(1 + N / a)
	}
	x = N / (N + a)
	return(c(a,x))
}

# fits the broken stick distribution
fitBS<-function(n)	{
	n = sort(n,decreasing=T)
	best = 999999
	spp = length(n)
	bestr = spp
	p = n / sum(n)
	for (i in spp:(spp * 3))	{
		q = rep(0,i)
		for (j in 1:i)	{
			for (k in 0:(i-j))	{
				q[j] = q[j] + 1 / ( i - k )
			}
			q[j] = q[j] / i
		}
		q = q / sum(q,na.rm=T)
		kl = sum(p * log(p / q[1:spp]))
		if (kl < best)	{
			best = kl
			bestr = i
			bestbs = q
		}
		if (bestr != i)	{
			break
		}
	}
	return(list(distribution=bestbs,richness=bestr,fit=best))
}

# fits the Zipf distribution
fitZipf<-function(n)	{
	n = sort(n,decreasing=T)
	spp = length(n)
	p = n / sum(n)
	co = lm(log(n) ~ log(1:spp),weights=sqrt(n))$coefficients
	z = co[2]
	best = 999999
	for (i in 1:100)	{
		lastexp = z
		z = z + rnorm(1,sd=0.01)
		q = exp(log(1:1000) * z)
		q = q / sum(q)
		kl = sum(p * log(p / q[1:spp]))
		if (kl < best)	{
			best = kl
			bestexp = z
			bestzipf = q[1:spp]
		}
		if (bestexp != z)	{
			z = lastexp
		}
	}
	return(list(distribution=bestzipf,exponent=as.numeric(bestexp),fit=best))
}

