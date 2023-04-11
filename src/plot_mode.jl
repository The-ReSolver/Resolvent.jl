# This file contains a number of plotting functions for channel flow modes.

# TODO: test these

function plot_spacetime_mode_at_time(spacetime_mode::AbstractVector{AbstractArray{Float64, 3}}, ti::Int)
    # extract grid of points
    grid = points(get_grid(spacetime_mode))

    # plot the desired temporal mode
    p = controuf(grid[2], grid[1], spacefreq_mode[1][:, :, ti])
    quiver!(p, grid[2], grid[1], spacefreq_mode[2][:, :, ti], spacefreq_mode[3][:, :, ti])

    # plot formatting
    plot!(p, xlabel="z", ylabel="y")

    return p
end

function plot_space_mode_at_time(spacetime_mode::AbstractVector{AbstractArray{Float64, 3}}, tis::Vector{Int})
    # extract grid of points
    grid = points(get_grid(spacetime_mode))

    # initialise vector to hold plots
    ps = Vector{Plots.Plot}(undef, length(tis))

    # plot the desired temporal mode
    for (i, ti) in enumerate(tis)
        ps[i] = controuf(grid[2], grid[1], spacefreq_mode[1][:, :, ti])
        quiver!(ps[i], grid[2], grid[1], spacefreq_mode[2][:, :, ti], spacefreq_mode[3][:, :, ti])

        # plot formatting
        plot!(ps, xlabel="z", ylabel="y")
    end

    return ps
end
