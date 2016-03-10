	

module funcapp




	# use chebyshev to interpolate this:
# f(x) = x + 2x^2 - exp(-x)$ for $x\in[-3,3]
	
	function q1(n)

n = 15	
using PyPlot
using FastGaussQuadrature: gausschebyshev
nodes = gausschebyshev(n)
nodesbis = 3nodes[1]

truef(x) = x + 2x.^2 - exp(-x) 

points = []
for i=1:n
	push!(points, truef(nodesbis[i]))
end
points = convert(Array{Float64,1},points)

phi = [cos((n-i+0.5)*(j-1)*pi/n) for i in 1:n,j in 1:n]
phi = convert(Array{Float64,2},phi)
coef = inv(phi)*points

## z = 2*((nodesbis + 3)/6) - 1
## T(y) = [cos(acos(y)*j) for j in 0:n-1]
# approx(y) = dot(coef, T(y/3))

function t(z::Real,j::Int)
	if j==0 return 1
	elseif j ==1 return z
	else return 2*z*t(z,j-1) - t(z,j-2)
	end 
end

approxrec(z) = sum([coef[i+1] *t(z/3,i) for i in 0:n-1])


v = linspace(-3,3, 50)
pointsbis = map((x)-> approxrec(x), v)

x = linspace(-3,3, 1000)
plot(x, truef(x), "b-")
plot(v, pointsbis, "r+")

m = 1000
w = linspace(-3,3, m)
approxrecvec(l) = map((s) -> approxrec(s), l)
truefvec(l) = map((s) -> truef(s),l)
resid(l) = abs(truefvec(l)-approxrecvec(l))
residmax = pop!(maximum(resid(w), 1))

if residmax > 1.0^(-9)
	println("Reject the approximation")
else println("Accept the approximation")
end


	end

	function q2(n)
		
using ApproxFun

S=Chebyshev([-3,3])
x=points(S,n)
v = x + 2x.^2 - exp(-x)
fun =Fun(ApproxFun.transform(S,v),S)
ApproxFun.plot(fun)

	end


	# plot the first 9 basis Chebyshev Polynomial Basisi Fnctions
	function q3()


x = linspace(-1,1,1000)

che = [cos(acos(x)j) for j in 0:8]

fig,axes = subplots(3,3,figsize=(10,5))

for i in 1:3
	for j in 1:3
        ax = axes[j,i]
        count = i+(j-1)*3
        ax[:plot](x,che[i+(j-1)*3])
        ax[:set_title]("Basis function $(count-1)")
        ax[:yaxis][:set_visible](false)
        ax[:xaxis][:set_visible](false)
        ax[:set_xlim](-1.0,1.0)
    end
end
fig[:canvas][:draw]()


	end

	ChebyT(x,deg) = cos(acos(x)*deg)
	unitmap(x,lb,ub) = 2.*(x.-lb)/(ub.-lb) - 1	#[a,b] -> [-1,1]

	type ChebyType
		f::Function # fuction to approximate 
		nodes::Union{Vector,LinSpace} # evaluation points
		basis::Matrix # basis evaluated at nodes
		coefs::Vector # estimated coefficients

		deg::Int 	# degree of chebypolynomial
		lb::Float64 # bounds
		ub::Float64

		# constructor
		function ChebyType(_nodes::Union{Vector,LinSpace},_deg,_lb,_ub,_f::Function)
			n = length(_nodes)
			y = _f(_nodes)
			_basis = Float64[ChebyT(unitmap(_nodes[i],_lb,_ub),j) for i=1:n,j=0:_deg]
			_coefs = _basis \ y  # type `?\` to find out more about the backslash operator. depending the args given, it performs a different operation
			# create a ChebyType with those values
			new(_f,_nodes,_basis,_coefs,_deg,_lb,_ub)
		end
	end
	
	# function to predict points using info stored in ChebyType
	function predict(Ch::ChebyType,x_new)

		true_new = Ch.f(x_new)
		basis_new = Float64[ChebyT(unitmap(x_new[i],Ch.lb,Ch.ub),j) for i=1:length(x_new),j=0:Ch.deg]
		basis_nodes = Float64[ChebyT(unitmap(Ch.nodes[i],Ch.lb,Ch.ub),j) for i=1:length(Ch.nodes),j=0:Ch.deg]
		preds = basis_new * Ch.coefs
		preds_nodes = basis_nodes * Ch.coefs

		return Dict("x"=> x_new,"truth"=>true_new, "preds"=>preds, "preds_nodes" => preds_nodes)
	end





	function q4a(deg=(5,9,15),lb=-5.0,ub=5.0)


# Define the function to approximate
fct(x) = 1.0./(1+25.*x.^2)
l = linspace(lb,ub,1000) # Points to extrapolate

fig,axes = subplots(1,2,figsize=(10,5))

colors   = ["blue","red","green"]

ax = axes[1,1]
ax[:set_title]("Chebyshev nodes")

for j in 1:3

	## You do not have to change the scales of the nodes, (there are already between -1 and 1)
	nodes = gausschebyshev(deg[j]+1)[1]*5
	C = ChebyType(nodes,deg[j],lb,ub,f)
	x = predict(C,l)["x"]
	y1 = predict(C,l)["truth"]
	y2 = predict(C,l)["preds"]

	ax[:plot](x,y1, color = "black")  # True function
	ax[:plot](x,y2, label="deg=$(deg[j])", color=colors[j])
	ax[:set_ylim](-0.1,1.1)
	ax[:legend](loc="upper right")	
end


ax = axes[2,1]
ax[:set_title]("Uniform nodes")

for j in 1:3

	nodes = linspace(-5,5,deg[j]+1)

	C = ChebyType(nodes,deg[j],lb,ub,f)
	x = predict(C,l)["x"]
	y1 = predict(C,l)["truth"]
	y2 = predict(C,l)["preds"]

ax[:plot](x,y1, color = "black") # True function
ax[:plot](x,y2, label="deg=$(deg[j])", color=colors[j])
ax[:set_ylim](-0.1,1.1)
ax[:legend](loc="upper right")	

end

fig[:canvas][:draw]()
println("With Chebyshev nodes, the appromixation improves as the number of nodes increases, but the overall fit is not that great with few points (especially in the middle)")
println("With uniform nodes, the approximation is a bit better to the Chebyshev with very few nodes but it gets worse at the boundary of the set when you increase the number of nodes")

	end

	function q4b()

using ApproXD
using Distributions
using PyPlot
runge(x) = 1.0./(1+25.*x.^2)
nknots = 13
deg = 3
lb=-5.0
ub=5.0
bs1 = BSpline(nknots,deg,lb,ub)

nevals = 5 * bs1.numKnots

		G(k,s) = GeneralizedPareto(k,s,0)
		pf(k,s) = quantile(GeneralizedPareto(k,s,0),linspace(0.05,cdf(G(0.5,1),5),6))
		myknots = vcat(-reverse(pf(0.5,1)),0.0,pf(0.5,1))

		bs2 = BSpline(myknots,deg)


		# get coefficients for each case

		eval_points = collect(linspace(lb,ub,nevals))  
		c1 = getBasis(eval_points,bs1) \ runge(eval_points)
		c2 = getBasis(eval_points,bs2) \ runge(eval_points)	# now \ implements the penrose-moore invers. see `pinv`. regression.

		# look at errors over entire interval
		test_points = collect(linspace(lb,ub,1000));
		truth = runge(test_points);
		e1 = getBasis(test_points,bs1) * c1 - truth;
		e2 = getBasis(test_points,bs2) * c2 - truth;
		figure(figsize=(8,7))
		subplot(211)
		plot(test_points,truth,lw=2)
		ylim(-0.2,1.2)
		grid()
		title("Runge's function")
		subplot(212)
		plot(test_points,e1,label="equidistant",color="blue")
		plot(test_points,e2,label="concentrated",color="red")
		plot(unique(bs1.knots),zeros(nknots),color="blue","+")
		plot(myknots,ones(nknots)*ylim()[1]/2,color="red","o")
		ylim(minimum(e1)-0.1,maximum(e1)+0.1)
		grid()
		legend(loc="upper right",prop=Dict("size"=>8))
		title("Errors in Runge's function")


### my answer for this question is not only from me, but I read it and understood it.


	end

	function q5()


###I had no time to cover this question
		
	end


		# function to run all questions
	function runall()
		println("running all questions of HW-funcapprox:")
		q1(15)
		q2(15)
		q3()
		q4a()
		q4b()
		q5()
	end


end

