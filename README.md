# gnss_multipath

Tutorial to understand the impact of multipath on GNSS signals

## Introduction

Multipath is a phenomenon that causes GNSS receiver synchronization errors. The objective of this tutorial is to understand its impact on the receiver. We will see that the magnitude of the multipath error depends on several factors such as:
- The amplitude, delay and phase of the multipath propagation with respect to the direct signal.
- GNSS receiver properties (chip spacing and front-end bandwidth)
- The modulation of the GNSS signal

## Multipath interference 

The multipath interference phenomenon occurs when the signal received at the antenna is the sum of :
- the direct signal, called line of sight (LOS)
- at least one echo due to the reflection/diffraction of the LOS on an obstacle.

The echo is characterized by 
- Its attenuation with respect to the LOS
- its code delay offset from the LOS
- its carrier phase shift with respect to the LOS.

## Non-line-of-sight (NLOS)

Non-line-of-sight (NLOS) occurs when:
- the received signal does not contain the direct signal (the satellite is masked by an obstacle)
- at least one echo due to the reflection/diffraction of the LOS on an obstacle (NLOS).

The propagation time measured is that of the echo, so the distance between the satellite and the antenna is overestimated.

NLOS can be detected by
- Monitoring the signal amplitude drop (C/N0)
- 3D map of the environment 
- Fisheye camera on the roof of the vehicle + attitude estimation sensors.
- Low altitude satellite rejection

## Multi-trip Envelope Assumptions

### Objective:

- Evaluate the maximum magnitude of error that can affect a receiver.
- Evaluate the robustness of a receiver or signal to multipath.

### Assumption:

- The received signal is the sum of a direct signal and a single echo.
- The echo is either in phase or in phase opposition with the direct signal.
- The amplitude of the echo is assumed to be constant
- The error due to multipath propagation is plotted as a function of the delay of the echo relative to the direct signal.

### Method:

- Define the autocorrelation function of the signal
- Add the correlation of the local replica and the echo
- Build the S-curve (discriminator)
- Find where the discriminator is equal to zero

### Disadvantages: 
- Does not take into account loop filtering effects
- Does not take into account the dynamics (evolution of the multipath error over time)