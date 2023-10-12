# TODO: add forces by walls

using StaticArrays, Random, Statistics, DataFrames
using BenchmarkTools
# using FLoops
using LinearAlgebra: norm
using CSV
# using Tables: table
using .Threads
using ThreadsX
using Distributed: @distributed, @sync, @everywhere
@everywhere using SharedArrays

using ProgressMeter


mutable struct Particle
  pos::MVector{2, Float64}
  vel::MVector{2, Float64}
  size::Float64
  damp::Float64
  # color::RGBAf
  density::Float64
end

struct Sim
  dt::Float64
  tend::Float64
  radius::Float64
  # x_domain#
  # y_domain
end


function update!(particles::Vector{Particle}, sim::Sim, condx::Float64, condy::Float64, pressureMultiplier::Float64, randomMultiplier::Float64)

  gravity = 00.0;

  # TODO: add parallelization to the SEPARATE "for" loops
  # TODO: add @sync @distributed i.e. 
  

  # calculate densities
  # apply gravity
  # @floop begin
  
  # ThreadsX.foreach(enumerate(particles)) do (particleIndex, p)
  # println("Update step")
  # for (particleIndex, p) in enumerate(particles)
  @sync @distributed for p in particles
    p.vel += [0.0,-1.0] * gravity * sim.dt;
    #particle_densities = [Density(particles, p.pos, smooth_radius) for p in particles]
    p.density = Density(particles, p.pos, sim.radius) 
  end
# end

# calculate pressure forces
# @floop begin
  # TODO: use verlet
  # TODO: add viscous forces (compute laplacian)
  forceArray = zeros(length(particles))
  for (particleIndex, p) in enumerate(particles)
    forces = Force(particles, particleIndex, sim.radius, pressureMultiplier)
    forces = forces / p.density
    p.vel = p.vel + forces * sim.dt; # no inertia (gas-like?)
    randVel = [rand()*rand([1.0,-1.0]), rand()*rand([1.0,-1.0])]
    p.vel += randVel / norm(randVel) * randomMultiplier
    
    forceArray[particleIndex] = norm(forces)
  end 
# end

# solve for collisions
# @floop for p in particles
# @floop begin 
  for p in particles

    p.pos += p.vel * p.damp * sim.dt;

    # check y
    if abs(p.pos[2]) > condy
      
      p.pos[2] = condy * sign(p.pos[2])
      p.vel[2] *= -1 * p.damp
      #scatter!(ax, transpose(p.pos), color=:blue)
    # check x
    end
    if abs(p.pos[1]) > condx
      p.pos[1] = condx * sign(p.pos[1])
      p.vel[1] *= -1 * p.damp
      #scatter!(ax, transpose(p.pos), color=:red )
    #else
    end
    # print("P($particleIndex) positions: $(p.pos)\t")
  end
  # println("")
  # println("-------------")
  # end
  return forceArray
  
end


#main()

function kernel(radius::Float64, distance::Float64)
  # later: implement conditional cluster
  # S. Clavet, pvfs.pdf
  if distance >= radius return 0 end;
  vol = π * radius^4 / 6
  #ρ = maximum([0, radius - dist*dist])
  W = (radius - distance)^2
  return W / vol
end

function kernelDerivative(radius::Float64, distance::Float64)
  # kernel derivative
  if distance >= radius return 0 end;
  scale = 12 / (radius^4 * π)
  #ρ = maximum([0, radius - dist*dist])
  #W = (radius - distance)^2
  return (distance - radius) * scale
end

@everywhere function Density(p::Vector{Particle},samplePoint::MVector{2,Float64}, radius::Float64)
  den = 0.0
  mass = 1.0
  
  for particle in p
    dst = norm((particle.pos - samplePoint))
    weight = kernel.(radius, dst)
    den += mass * weight
  end
  return den
end

# function A_func(p::Vector{Particle}, A_j::Vector,ρ_j::Vector, sample_point::Vector{Float64})
#   # calculates field operator (A) over sampled particle points 
#   # with (A_j) from Koschier's paper
#   # ρ_j: density vector
  
#   A = 0.0
#   mass = 1.0
#   for particle in p
#     dst = norm((particle.pos - sample_point))
#     weight = kernel.(radius, dst)
#     A += A_j * mass / ρ_j * weight
#   end
#   return A

# end

function randomDir()
  return rand(2)
end


function densityToPressure(ρ::Float64, pressureMultiplier::Float64)
  targetDensity = 1.0
  # pressureMultiplier = 50.0

  densityError = ρ - targetDensity
  pressure = densityError * pressureMultiplier
  return pressure
  
end

function Force(p::Vector{Particle}, particleIndex::Int, radius::Float64, pressureMultiplier::Float64)
  # calculates force only between particles
  # have to make sure no singularity is present in the distances
  # we do this by
    # Iterating over i≠j
    # Inserting a Lennard-Jones potential (later)
      # for now, just add a random velocity vector

  mass = 1.0
  force = zeros(2)
  particleDensity = p[particleIndex].density
  particlePressure = densityToPressure(particleDensity, pressureMultiplier)
  # println("")
  # println("Force Calc")
  for (otherParticleIndex, particle) in enumerate(p)
    # sample point = particle i
    # rest of points = particle j, j≠i
    
    # otherDensity = densityVec[otherParticleIndex]
    otherDensity = particle.density;

    if particleIndex==otherParticleIndex 
      #println("id: $otherParticleIndex")
      continue; 
    end;
    direction = p[particleIndex].pos - p[otherParticleIndex].pos
    # println("id: $otherParticleIndex, direction: $direction")
    distance  = norm(direction)
    
    # direction = dir/dst
    #
    # check if distance is zero
    if distance == 0

      direction = randomDir()
    else
      #println(distance)
      direction = direction / distance
    end

    dWdx = kernelDerivative(radius, distance);
    otherPressure = densityToPressure(otherDensity, pressureMultiplier)
    # calculate shared force
    force += -(particlePressure + otherPressure) / 2 * direction * dWdx * mass / otherDensity;
    # force += otherPressure * direction * dWdx * mass / otherDensity;
  end
  return force
end


#begin
#begin #main()
function main(n_particles::Int, smooth_radius::Float64, pressure_multiplier::Float64, dampen_factor::Float64)
  # convention: snake_case for variables
  #             camelCase for passed values in functions
  RNG = MersenneTwister(1223)
  x0 = 1.0
  y0 = 1.0
  dx = 0.2
  n_dom = 100
  #n_particles=500
  # dampen_factor= dampen_factor / 
  particle_size = 0.001
  
  #smooth_radius = 0.5
  T = Sim(0.01, 10.0, smooth_radius)
  
  bounds_size = [2*x0, 2*y0]
  half_bound_size = bounds_size / 2 - ones(2) * particle_size;




  # x_domain = vcat(collect(-n_dom/2:0.0)*2/n_dom * (dx) .- (x0 - dx),
  #                 collect(0.0:n_dom/2)*2/n_dom  * (dx) .+ (x0 - dx))
  # x_domain = (-n_dom/2:n_dom/2) * 2/n_dom * x0/3
  # y_domain = (-n_dom/2:n_dom/2) * 2/n_dom * y0/3
  x_domain = (-n_dom:0) * 1/n_dom * dx .- (x0 - dx) .+ dx; y_domain = (-n_dom:0) * 1/n_dom * dx .- (y0 - dx) .+ dx
  # x_domain = (-n_dom/2:n_dom/2) * 2/n_dom * x0; y_domain = (-n_dom/2:n_dom/2) * 2/n_dom * y0

  r_dom() = MVector{2}([rand(x_domain),rand(y_domain)])
  v_dom() = MVector{2}([rand(), rand()])
  particles = @showprogress [Particle(r_dom(), v_dom(),particle_size,dampen_factor, 0.0) for _ in 1:n_particles]
 
  # all samples
  # fluid = [Density(particles, [x,y], T.radius) for x in x_domain, y in y_domain]
  # TODO: calculate fluid in real time for the animation instead
  
  
  # isfile(cloud_file) ? rm(file_path) : 0;
  # dfx = DataFrame([name => Vector{Float64}(undef, 2) for name in ["x"]]) # empty dataframe
  t_range = 0:T.dt:T.tend
  dfx = zeros(length(t_range), n_particles)
  dfy = zeros(length(t_range), n_particles)
  df_force = zeros(length(t_range), n_particles)
  
  # run simulation
  # ThreadsX.foreach(enumerate(t_range)) do (row, t)
  random_multiplier=2.0
  @showprogress for (row, t) in enumerate(t_range)
    forces = update!(particles, T, half_bound_size[1], half_bound_size[2], pressure_multiplier, random_multiplier)
    positions = transpose(mapreduce((p::Particle) -> p.pos, hcat,  particles)) # get positions

    dfx[row,:] = @view positions[:,1]
    dfy[row,:] = @view positions[:,2]
    df_force[row,:] = forces
    #println("t=$t")
    #densities = transpose(mapreduce((p) -> p.density, hcat,  particles))
    #fluid = [Density(particles, [x,y], T.radius) for x in x_domain, y in y_domain]
    #push!
    #dfx = DataFrame(t=positions[:,1])
    #dfy = DataFrame(y=positions[:,2],t=t)
    # CSV.write(file_path_x, dfx, append=true)
    # CSV.write(file_path_y, dfy, append=true,newline=';', header=false, transpose=true)
    # CSV.write(cloud_path, df, append=true,newline='\n', header=false)
  end

  names = vcat(["n$i" for i in 1:n_particles], "t")
  dframe_x = DataFrame([dfx t_range], names)
  dframe_y = DataFrame([dfy t_range], names)
  dframe_f = DataFrame(df_force, names[1:end-1])
  # # create for storage  
  file_path_x = "./x_pos.csv"
  file_path_y = "./y_pos.csv"
  file_path_f = "./forces.csv"
  isfile(file_path_x) ? rm(file_path_x) : 0;
  isfile(file_path_y) ? rm(file_path_y) : 0;
  isfile(file_path_f) ? rm(file_path_f) : 0;
  
  CSV.write(file_path_x, dframe_x)
  CSV.write(file_path_y, dframe_y)
  CSV.write(file_path_f, dframe_f)

  # save environment variables

  #CSV.write("env.csv", )

end
# let
  
#   fig = Figure(resolution = (1200,600))
#   ax = Axis(fig[1,1])
#   #scatter!(ax, transpose(sample))
#   hm=heatmap!(ax, x_domain, y_domain, fluid, colormap=:bam)
#   Colorbar(fig[1,2],hm)

#   positions = transpose(mapreduce((p) -> p.pos, hcat,  particles)) # get positions
#   force_mat = mapreduce(permutedims, vcat, forces)  
  
#   scatter!(ax, positions, color=:black, markersize=10)
#   arrows!(ax, positions[:,1],positions[:,2], force_mat[:,1], force_mat[:,2], lengthscale=0.004)
#   #lines!(ax,smooth_radius*cos.(θ) .+ sample[1], smooth_radius*sin.(θ) .+ sample[2])
#   ax.aspect=DataAspect()
#   #colsize!(fig.layout, 1, Aspect(1, 1))    
#   #println(sum(fluid[:]))
#   fig
# end


# # animation
# begin
#   X = CSV.read("./x_pos.csv", DataFrame, header=true)
#   Y = CSV.read("./y_pos.csv", DataFrame, header=true)
  
#   fig = Figure()
#   ax = Axis(fig[1,1])
#   [scatter!(ax, X[:,p], Y[:,p],color=(:black, 0.5), markersize=5) for p in 1:size(X,2)-1]
#   fig
# end
