# This file contains the definitions required for a general mode implementation
# (probably to moved to the interface package at a later date) which defines
# the basis to express flows upon.

#=
* A mode is simply a function onto which a flow field can be projected onto.

* The goal of any mode type is to provide this basic functionality, not only
* for a single mode, but for collection of modes forming a basis (a simple
* vector would do for this for any flow type).

* This can be expressed abstractly, thus allowing a consistent interface to be
* defined for all the modal analysis techniques used.

* It will also be useful to be able to project modes onto each other, mainly
* for the sake of comparison of the modes generated from different techniques.

* Therefore, the following needs to be implemented in this file:
*   - an abstract mode type
*   - projection methods (for both flow fields and other modes) on this abstract type
*   - a similar projection method for a collection of modes (inherits from the abstract method so shouldn't need any concrete implementations)
*   - a concrete implementation of this for channel flow (restricts us to modes as a wall-normal profile)
*   - a concrete implementation of the projection method
*   - plotting method for the modes

* Are there other important operations that can be expressed in terms of a few simple abstract operations???
=#


