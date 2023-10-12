using CairoMakie, CSV, DataFrames
using LinearAlgebra: norm
using ProgressMeter

struct Sim
  dt::Float64
  tend::Float64
  radius::Float64
  # x_domain#
  # y_domain
end

function kernel(radius::Float64, distance::Float64)
  # later: implement conditional cluster
  # S. Clavet, pvfs.pdf
  if distance >= radius return 0 end;
  vol = π * radius^4 / 6
  #ρ = maximum([0, radius - dist*dist])
  W = (radius - distance)^2
  return W / vol
end

function Density(X::Vector{Float64},Y::Vector{Float64},samplePoint::Vector{Float64}, radius::Float64)
  den = 0.0
  mass = 1.0
  
  #ThreadsX.foreach(zip(X, Y)) do (x,y)
  for (x,y) in zip(X,Y)
    dst = norm([x,y] - samplePoint)
    weight = kernel.(radius, dst)
    den += mass * weight
  end
  return den
end

#function Density(X::Observable, Y::Observable)


function animate(percentage)

  begin 
    xdata = CSV.read("./x_pos.csv", DataFrame); 
    ydata = CSV.read("./y_pos.csv", DataFrame) 
  end;

  x0 = 2.0
  y0 = 1.0
  dr = 0.01
  radius = 0.5
  particle_size = 0.001

  xdom = Observable(-x0:dr:x0)
  ydom = Observable(-y0:dr:y0)
  
  
  #row_arr = (1:20)
  # FIGURE ##########################
  fig =  Figure(theme=theme_black())
  ax = Axis(fig[1,1], xlabel="x", ylabel="y")
  X = Observable(collect(xdata[1,begin:end-1]))
  Y = Observable(collect(ydata[1,begin:end-1]))
  
  f(i,j) = Density(X.val,Y.val, [i,j], radius)
  
  fluid = f.(xdom.val, ydom.val')
  fluid_obs = Observable(f.(xdom.val, ydom.val'))
  joint_limits = (minimum(fluid),maximum(fluid))
  
  hm=heatmap!(ax, xdom, ydom, fluid_obs, colormap=:curl,colorrange = joint_limits)
  scatter!(ax, X, Y, color=:black)
  ax.aspect=DataAspect()
  # cb=Colorbar(fig[1,2], hm, colorrange = (min_col, max_col))
  #  resize!(f.scene, gridlayoutsize(f.layout) .+ (32, 32))
  resize_to_layout!(fig)
  
  
  framerate = 50
  all_rows = convert(Int, round(size(xdata, 1) * percentage/100))
  prog = Progress(all_rows)
  record(fig, "time_animation.mp4", 1:all_rows; framerate = framerate) do row
    X[] = collect(xdata[row,begin:end-1])
    Y[] = collect(ydata[row,begin:end-1])
    fluid_obs[] = f.(xdom.val, ydom.val')
    next!(prog)
  end
end

function compareStates()
  begin 
    xdata = CSV.read("./x_pos.csv", DataFrame); 
    ydata = CSV.read("./y_pos.csv", DataFrame) 
  end;
  x0 = 2.0
  y0 = 1.0
  dr = 0.01
  radius = 0.5
  particle_size = 0.001

  xdom = (-x0:dr:x0)
  ydom = (-y0:dr:y0)
  
  
  #row_arr = (1:20)
  # FIGURE ##########################
  fig =  Figure(theme=theme_black())
  ax = Axis(fig[1,1], xlabel="x", ylabel="y")
  X = (collect(xdata[1,begin:end-1]))
  Y = (collect(ydata[1,begin:end-1]))
  
  f(i,j) = Density(X,Y, [i,j], radius)
  
  fluid = f.(xdom, ydom')
  fluid_obs = (f.(xdom, ydom'))
  joint_limits = (minimum(fluid),maximum(fluid))
  
  hm=heatmap!(ax, xdom, ydom, fluid_obs, colormap=:curl,colorrange = joint_limits)
  scatter!(ax, X, Y, color=:black)
  ax.aspect=DataAspect()
  # cb=Colorbar(fig[1,2], hm, colorrange = (min_col, max_col))
  #  resize!(f.scene, gridlayoutsize(f.layout) .+ (32, 32))
  resize_to_layout!(fig)
  

  fig2 =  Figure(theme=theme_black())
  ax2 = Axis(fig2[1,1], xlabel="x", ylabel="y")
  X = (collect(xdata[end,begin:end-1]))
  Y = (collect(ydata[end,begin:end-1]))
  
  f(i,j) = Density(X,Y, [i,j], radius)
  
  fluid = f.(xdom, ydom')
  fluid_obs = (f.(xdom, ydom'))
  
  hm=heatmap!(ax2, xdom, ydom, fluid_obs, colormap=:curl,colorrange = joint_limits)
  scatter!(ax2, X, Y, color=:black)
  ax2.aspect=DataAspect()
  # cb=Colorbar(fig[1,2], hm, colorrange = (min_col, max_col))
  #  resize!(f.scene, gridlayoutsize(f.layout) .+ (32, 32))
  resize_to_layout!(fig2)
  return fig, fig2
end