---
layout: math_post
title: Fourth order Runga-Kutta method
author: Jaco Stroebel
permalink: fourth-order-runga-kutta
tags: [python, runga-kutta, initial value problems]
---

When modelling phenomena we often want a way of describing change related to the phenomena, like population growth, particle speed, etcetera.
More often than not these differential equations <?A differential equation is any function that is a derivative of another function, for example describing the rate of change[http://en.wikipedia.org/wiki/Differential_equation]> are what is called non-elementary, meaning that they cannot be solved using trivial (analytical?) integration methods.
For such problems we have to use computational approximation.
There are various algorithms available for these approximations, some of them are trivial to reproduce by hand, some not so much.

####Theory

Let us consider the initial value problem

$$
\frac{\text{d}y}{\text{dt}} = f(\text{t}, y(\text{t})) (should we use \Delta?)
\; \text{with} \; y(\tau) = b
$$

The function `\(y(\text{t})\)` is the unknown function that we would like to approximate, we can extrapolate the derivative `\(\frac{\text{d}y}{\text{dt}}\)` from observation.
Approximation with the Runga-kutta method cannot give us a closed from of the function, only the value at certain points.
It is said to be a discrete approximation.

Let us approximate the value for `\(y(\text{T})\)`.
Divide the interval
`\([\tau, \text{T}]\)`
into **N** subintervals with length
`\(\text{h} = \frac{\text{T} - \tau}{\text{N}}\)`.

$$
\tau = \text{t}\_0 \lt \text{t}\_1 \lt \cdots \lt \text{t}\_N = \text{T}
$$

In the subinterval `\([t_{k+1}, t_k]\)` it follows that we can express the difference in the values of the function `\(y(\text{t})\)` at the boundries like so, 

$$
y(\text{t}\_{k+1}) - y(\text{t}\_k) = \int\_{t\_k}^{t\_{k+1}} f(\text{t}, y(\text{t})) \text{dt}
$$ 

It was already mentioned that we cannot integrate the function, so we make use of a numerical integration technique.
The technique is called the 1/3 Simpson's Rule, credited to the English mathematician Thomas Simpson (1710-1761).
With this technique we find three points in the interval `\([t_{k+1}, t_k]\)`, `\(t_k\)`,`\(\frac{t_{k+1} - t_k}{2}\)` and `\(t_{k+1}\)`.
From this we can construct a Langranian polynomial of degree 2 that we can easily integrate.
Now we have 

$$
y(\text{t}\_{k+1}) - y(\text{t}\_k) \approx
$$
`\(
\frac{h}{6}\left(f\left(\text{t}_k, y\left(\text{t}_k\right)\right) + 4f\left(\text{t}_k + \frac{h}{2}, y\left(\text{t}_k + \frac{h}{2}\right)\right) + f\left(\text{t}_{k+1}, y\left(\text{t}_{k+1}\right)\right) \right)
\)`

From here we add `\(y(\text{t})\)` on both sides and substitute in the appropriate function values then we have the approximation that we want.
The snag however, is that we do not have the function values since we do not know what the function is.
We solve this problem by approximating those values and using each successive predictive value to approximate the next.
The shorter we take our steps, that is the larger **N** is, the more accurate our approximation will be. 

####Implementation 

Let `\(y'(\text{t}) = \text{t}y + y + \text{t}^2\)` <? this is the differential equation (function?) that we want to approximate> with `\(y(0) = 2\)` <? a known value at point 0>.
We can approximate any arbitrary `\(\text{T}\)`, in this case it will be 0.4.
(any arbitrary y(T), in this case for T=0.4)?


```Python
def f(y, t):
    return t*y + y + t**2

tau = 0 # initial t
T = 0.4 # final t
N = 100 # number of steps
h = (T - tau) / N # step size

t = [0] # t in each step
y = [2] # results

for k in range(0, N):
    F_1 = f(y[k], t[k]) <?call our mysterious function with the first (or current) result and current t>
    F_2 = f((y[k] + h/2*F_1), (t[k] + h/2))
    F_3 = f((y[k] + h/2*F_2), (t[k] + h/2))
    F_4 = f((y[k] + h*F_3), (t[k] + h))
    
    y.append(y[k] + h/6*(F_1 + 2*F_2 + 2*F_3 + F_4))
    <? in other words, the result consists of 6 parts, one part F_1,
    two parts each F_2 and F_3 and one part F_4, all scaled by the step size>
    
    t.append(t[k] + h)

print(y[N])
```

The variables `F_1`, `F_2`, `F_3` and `F_4` are what we refer to as predictive values, the values approximate the values for `\(y\)` that we cannot explicitly find.
Each of these values increase in distance from our initial value.
The longer an approximation step is, the larger the error.
That is why we use the previous approximation to approximate the next.

The loop iterates over the amount of steps that we would like to perform, each time appending the approximate value of `\(y\)` to the list `y`. When it is done and we want to find the value of `\(y(\text{T})\)`, we print the value of `y[N]`. 
That is because there are N + 1 elements in `y`, so it is a perfectly valid access operation.

This spesific implementation should complete in linear time.

####Done

Fourth Order Runga-Kutta is a very accurate method for the numerical approximation of ODE's, it is also fairly simple to implement. Expanding this process to systems of equations is as simple as defining another function and then approximating the predictive values of both functions in the same loop.
For higher derivatives can be expressed as a system of derivatives and can therefore also be computed using this method.
