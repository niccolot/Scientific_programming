# **2 body problem with central force**

This script is a numerical integration algorithm using $2^{\circ}$ order Runge-Kutta of a system composed of 2 bodies interacting with a coulomb-like force $\propto \displaystyle\frac{1}{r^2}$ both inside a central force field, like two planets orbiting around a fixed sun and interacting with themselves.

$$
\begin{cases}
\ddot{\vec{r_a}}(t) = -\displaystyle\frac{\vec{r_a}(t)}{r_a^3(t)} + \mu_b\displaystyle\frac{\vec{r_{ab}}(t)}{r_{ab}^3(t)}\\
\ \\
\ddot{\vec{r_b}}(t) = -\displaystyle\frac{\vec{r_b}(t)}{r_b^3(t)} - \mu_b\displaystyle\frac{\vec{r_{ab}}(t)}{r_{ab}^3(t)}
\end{cases}
$$

with $\vec{r_{ab}}(t) \equiv \vec{r}_b(t) - \vec{r}_a(t)$ the vector going from A to B. 

If one considers a planetary system this is equal of setting the sun's mass and coupling constant $M, G = 1$ and considering the reduced masses $\mu_a \equiv m_a/m_b$, $\mu_b \equiv m_b/m_a$