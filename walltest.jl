using CairoMakie, StaticArrays, Random, Statistics



function initialize(N::Int, xlims::SArray, ylims::SArray)
  #r = [[rand(RNG, xlims[1]:0.01:xlims[2]), rand(RNG, ylims[1]:0.01:ylims[2])] for _ in 1:N]
  #v = [[rand(RNG), rand(RNG)] for _ in 1:N]
  rx = rand(RNG, xlims[1]:0.01:xlims[2], N)
  ry = rand(RNG, ylims[1]:0.01:ylims[2], N)
  v = rand(RNG, N,2)
  return  [rx ry], v
end

function frame()
  # return frame limits
  ymin = -1.0; ymax = 1.0;
  xmin = -1.0; xmax = 1.0;
  return (xmin, xmax), (ymin, ymax)
end

begin
  RNG = MersenneTwister(123)
  n_particles = 3;
  particle_size = 0.1;
  dampen_factor = 1.0 - 0.1;
  boundsSize = 1.0
  halfBoundsSize = boundsSize / 2 .- ones(2) * particle_size;

  gravity = -1.0
  x_range = SVector(-1.0, 1.0); 
  y_range = SVector(-1.0, 1.0); 
  
  r, v = initialize(n_particles,x_range, y_range)

  xlim, ylim = frame()

  dt = 0.1
  T_end = 1.5;
  t_vec = 0:dt:T_end;
  steps = length(t_vec)
  
  r_vec = zeros(size(r,1),size(r,2), steps)
  v_vec = zeros(size(v,1),size(v,2), steps)
  
  v_down = [zeros(n_particles) ones(n_particles)]
  
  r_vec[:,:,1] = r
  
  #v_vec[:,:,1] = v
  
  for t in eachindex(t_vec)[2:end]
    v_vec[:,:,t] = v_vec[:,:,t-1] + v_down * gravity * dt 
    r_vec[:,:,t] = r_vec[:,:,t-1] + v_vec[:,:,t-1] * dt;
    # resolve collision
    #if abs(r_vec[])

  end  

  idr = CartesianIndices(r_vec)
  fig = Figure()
  axs = Axis(fig[1,1])
  [scatter!(axs, r_vec[idr[p,:,:]]) for p in 1:n_particles]
  fig
end