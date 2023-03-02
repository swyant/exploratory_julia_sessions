### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 4b2869dc-b92f-11ed-1059-eda4bc903236
begin 
	mutable struct testStruct
		f1::Int64
		f2::Int64
	end
end

# ╔═╡ 16ac4927-2032-4c69-a76b-3e89e2009342
begin
	function modStruct!(x::testStruct)
		x.f1 = x.f1 +2
		x.f2 = x.f2 +2 
		return
	end
end

# ╔═╡ 8830db69-2906-4052-85d3-75835bd3033f
baz = testStruct(2,2)

# ╔═╡ 99ea0b03-1999-43e2-b269-3a6b30ff602f
baz.f2 = 3

# ╔═╡ e7ace81e-5274-4c58-801b-bbac85d16af4
modStruct!(baz)

# ╔═╡ 43e0a2e1-b3ca-4355-859c-6ea27839e671
md"Very clear that the above two cells do not trigger re-evaluation of the cell below"

# ╔═╡ 6e0d838c-db01-48c5-97e7-0ffc96f590b0
baz.f2

# ╔═╡ Cell order:
# ╠═4b2869dc-b92f-11ed-1059-eda4bc903236
# ╠═16ac4927-2032-4c69-a76b-3e89e2009342
# ╠═8830db69-2906-4052-85d3-75835bd3033f
# ╠═99ea0b03-1999-43e2-b269-3a6b30ff602f
# ╠═e7ace81e-5274-4c58-801b-bbac85d16af4
# ╠═43e0a2e1-b3ca-4355-859c-6ea27839e671
# ╠═6e0d838c-db01-48c5-97e7-0ffc96f590b0
