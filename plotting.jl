using CairoMakie, CSV, DataFrames
using LinearAlgebra: norm
using ProgressMeter
using .Threads

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
  @inbounds for (x,y) in zip(X,Y)
    @fastmath dst = norm([x,y] - samplePoint)
    weight = kernel.(radius, dst)
    den += mass * weight
  end
  return den
end

#function Density(X::Observable, Y::Observable)


function animate(percentage, cmap, radius::Float64, dr::Float64)

  begin 
    xdata = CSV.read("./x_pos.csv", DataFrame); 
    ydata = CSV.read("./y_pos.csv", DataFrame) 
    fdata = CSV.read("./forces.csv", DataFrame)
  end;
  fdata = Matrix(fdata)

  x0 = 1.0
  y0 = 1.0
  #dr = 0.05
  # radius = 0.5
  particle_size = 0.001

  xdom = Observable(-x0:dr:x0)
  ydom = Observable(-y0:dr:y0)
  
  all_rows = convert(Int, round(size(xdata, 1) * percentage/100))
  fluid_vec = zeros(length(xdom.val), length(ydom.val), all_rows)
  
  f(x,y,i,j) = Density(x,y, [i,j], radius)

  @time XV = [collect(xdata[row,1:end-1]) for row in 1:all_rows]
  @time YV = [collect(ydata[row,1:end-1]) for row in 1:all_rows]
  # @time map(f, Mat[])
  prog = Progress(all_rows)

  @time begin
    @threads for row in 1:all_rows # this takes a lot of time
      @inbounds @fastmath fluid_vec[:,:,row] = f.(Ref(XV[row]),Ref(YV[row]), xdom.val, ydom.val')
      next!(prog)
    end
  end
  #row_arr = (1:20)
  # FIGURE ##########################
  fig =  Figure(theme=theme_dark())
  ax = Axis(fig[1,1], xlabel="x", ylabel="y", limits=(-x0, x0, -y0, y0))
  X = Observable(collect(xdata[1,begin:end-1]))
  Y = Observable(collect(ydata[1,begin:end-1]))
  
  # f(i,j) = Density(X.val,Y.val, [i,j], radius)
  
  fluid_obs = Observable(fluid_vec[:,:,1])
  #joint_limits = (minimum(fluid_obs.val),maximum(fluid_obs.val))
  
  hm=heatmap!(ax, xdom, ydom, fluid_obs, colormap=:curl,interpolate=true)
  C = Observable(fdata[1,:]);
  
  scatter!(ax, X, Y, color=C, colormap=cmap)
  ax.aspect=DataAspect()
  #cb=Colorbar(fig[1,2], hm)
  # resize!(f.scene, gridlayoutsize(f.layout) .+ (32, 32))
  resize_to_layout!(fig)
  
  
  framerate = 60
  prog = Progress(all_rows)
  anim_name = "time_animation.mp4"
  record(fig, anim_name, 1:all_rows; framerate = framerate) do row
    X[] = collect(xdata[row,begin:end-1])
    Y[] = collect(ydata[row,begin:end-1])
    C[] = fdata[row,:]
    fluid_obs[] = fluid_vec[:,:,row]
    next!(prog)
  end
  run(`mpv --loop-file=yes time_animation.mp4 `)
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