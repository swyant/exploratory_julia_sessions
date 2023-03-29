### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 41c2564a-ce6b-11ed-0ac3-2752d1deec22
using AtomsIO

# ╔═╡ bf7a1e75-61f9-4255-9668-6f3a25726f1a
using AtomsIOPython

# ╔═╡ 7f0b1f32-3cac-4707-b521-b770b451c74f
using AtomsBase

# ╔═╡ 985dfa48-53e2-4573-bddb-fac174e2b36b
using Unitful

# ╔═╡ 0b718e36-3aea-4359-8414-26416277893a
md"Raw Approach"

# ╔═╡ 8279854f-cc71-4f95-9d22-107085ac9886
begin 
	box1 = [[21.105598636581512,0.0,0.0],
	        [0.0,21.105598636581512,0.0],
	        [0.0,0.0,126.63359181948998]]u"Å"
end

# ╔═╡ 8d151a24-3566-4eb3-af55-0d56a5701d3a
bcs1 = [Periodic(), Periodic(), Periodic()]

# ╔═╡ aac6c11a-7e01-48db-8042-74fd13316ae9
lines = readlines("../ref_config/raw_coords")

# ╔═╡ 4e50d8d2-2dcb-4c0c-a4cb-93c33d89998f
pos1  = [parse.(Float64, split(li)) for li in lines]

# ╔═╡ 6d035aa0-48b3-4965-bbb4-7c810eedf834
atoms1 = [Atom(:Ar, ri * u"Å") for ri in pos1]

# ╔═╡ b60a0a84-b2ca-4f9e-b389-352483b1f65c
sys1 = FlexibleSystem(atoms1,box1,bcs1)

# ╔═╡ f514edad-b5af-4a0b-99ec-d0af6b8b7423
md"Reading standard xyz file -- Appears to fail"

# ╔═╡ fea940f8-51a0-4069-9742-4e4269bf853f
#system2 = load_system("../ref_config/dump_final.xyz")

# ╔═╡ 38abf557-84c9-4d6c-b7e4-a039d08a5309
#system2 = load_system("../ref_config/fixed_dump_final.xyz")

# ╔═╡ 0d9e0de7-de2c-4040-ae2b-5adfbcd45987
md"Reading (manually created) extxyz file"

# ╔═╡ f7c445db-6a5c-43e6-bce0-6ea99cf11b37
system3 = load_system("../ref_config/dump_final.extxyz")

# ╔═╡ b85218e6-5d9d-48e4-a36f-88b02404bed7
md"Reading (manually created) extxyz file with forces"

# ╔═╡ dca85534-2f5b-4ab3-b517-2853f797faaa
system4 = load_system("../ref_config/dump_final_wforces.extxyz")

# ╔═╡ 5bd3263f-85a2-4fae-b3c8-62ebd5d01736
atomkeys(system4)

# ╔═╡ d9c11b88-cd94-4c4d-b236-e8dad58761fe
system4[1][:forces]

# ╔═╡ 7b78e10f-7928-483f-ba4b-1f7151d6a360
md"try reading in the .lammps file via ASE"

# ╔═╡ 458c33b0-2ddc-4b2c-84a4-d5c898ad2bd8
system5 = load_system("../ref_config/dump_final_forces.lammps")

# ╔═╡ 0f52ae1c-c683-4aaa-a1db-a570bb2b6a8f
system5[1]

# ╔═╡ 0977492e-c906-493b-8370-06505f2f5692
system6 = load_system("../ref_config/fixed_dump_final_forces.lammps")

# ╔═╡ 564f902d-2d0c-4b80-81eb-03ace2cc729f
atomkeys(system6)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AtomsBase = "a963bdd2-2df7-4f54-a1ee-49d51e6be12a"
AtomsIO = "1692102d-eeb4-4df9-807b-c9517f998d44"
AtomsIOPython = "9e4c859b-2281-48ef-8059-f50fe53c37b0"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
AtomsBase = "~0.3.2"
AtomsIO = "~0.2.0"
AtomsIOPython = "~0.1.0"
Unitful = "~1.12.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "f22d3bc40d98eb07fb44d6b1b83bd8ec930cd060"

[[deps.ASEconvert]]
deps = ["AtomsBase", "CondaPkg", "PeriodicTable", "PythonCall", "Unitful", "UnitfulAtomic"]
git-tree-sha1 = "9acd5754c9ffe7acd932b377d239a07984806e74"
uuid = "3da9722f-58c2-4165-81be-b4d7253e8fd2"
version = "0.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AtomsBase]]
deps = ["PeriodicTable", "Printf", "StaticArrays", "Unitful", "UnitfulAtomic"]
git-tree-sha1 = "570c83c85ad2580ff0dc317da72766d8e4a69a70"
uuid = "a963bdd2-2df7-4f54-a1ee-49d51e6be12a"
version = "0.3.2"

[[deps.AtomsIO]]
deps = ["AtomsBase", "Chemfiles", "ExtXYZ", "Logging", "PeriodicTable", "Reexport", "Unitful", "UnitfulAtomic"]
git-tree-sha1 = "f3c8c4260903a99c798461915a68683c493bbfe2"
uuid = "1692102d-eeb4-4df9-807b-c9517f998d44"
version = "0.2.0"

[[deps.AtomsIOPython]]
deps = ["ASEconvert", "AtomsIO", "Reexport"]
git-tree-sha1 = "c5b470fadfc5e572be2af3b36a475b1cebdd7977"
uuid = "9e4c859b-2281-48ef-8059-f50fe53c37b0"
version = "0.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Chemfiles]]
deps = ["Chemfiles_jll", "DocStringExtensions"]
git-tree-sha1 = "9126d0271c337ca5ed02ba92f2dec087c4260d4a"
uuid = "46823bd8-5fb3-5f92-9aa0-96921f3dd015"
version = "0.10.31"

[[deps.Chemfiles_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d4e54b053fc584e7a0f37e9d3a5c4500927b343a"
uuid = "78a364fa-1a3c-552a-b4bb-8fa0f9c1fcca"
version = "0.10.3+0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.CondaPkg]]
deps = ["JSON3", "Markdown", "MicroMamba", "Pidfile", "Pkg", "TOML"]
git-tree-sha1 = "741146cf2ced5859faae76a84b541aa9af1a78bb"
uuid = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
version = "0.2.18"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "89a9db8d28102b094992472d333674bd1a83ce2a"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.ExtXYZ]]
deps = ["AtomsBase", "PeriodicTable", "StaticArrays", "Unitful", "UnitfulAtomic", "extxyz_jll"]
git-tree-sha1 = "86f941b3596a8842852e8a2a9b7482f2079accf6"
uuid = "352459e4-ddd7-4360-8937-99dcb397b478"
version = "0.1.11"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "SnoopPrecompile", "StructTypes", "UUIDs"]
git-tree-sha1 = "84b10656a41ef564c39d2d477d7236966d2b5683"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.12.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.MicroMamba]]
deps = ["Pkg", "Scratch", "micromamba_jll"]
git-tree-sha1 = "a6a4771aba1dc8942bc0f44ff9f8ee0f893ef888"
uuid = "0b3b1443-0f03-428d-bdfb-f27f9c1191ea"
version = "0.1.12"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d78db6df34313deaca15c5c0b9ff562c704fe1ab"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.5.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.PeriodicTable]]
deps = ["Base64", "Test", "Unitful"]
git-tree-sha1 = "5ed1e2691eb13b6e955aff1b7eec0b2401df208c"
uuid = "7b2266bf-644c-5ea3-82d8-af4bbd25a884"
version = "1.1.3"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PythonCall]]
deps = ["CondaPkg", "Dates", "Libdl", "MacroTools", "Markdown", "Pkg", "REPL", "Requires", "Serialization", "Tables", "UnsafePointers"]
git-tree-sha1 = "f27dabb05ec811675a9eefe49325a14ae7266b0b"
uuid = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
version = "0.9.12"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "b8d897fe7fa688e93aef573711cb207c08c9e11e"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.19"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "ca4bccb03acf9faaf4137a9abc1881ed1841aa70"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.10.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "bb37ed24f338bc59b83e3fc9f32dd388e5396c53"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.12.4"

[[deps.UnitfulAtomic]]
deps = ["Unitful"]
git-tree-sha1 = "903be579194534af1c4b4778d1ace676ca042238"
uuid = "a7773ee8-282e-5fa2-be4e-bd808c38a91a"
version = "1.0.0"

[[deps.UnsafePointers]]
git-tree-sha1 = "c81331b3b2e60a982be57c046ec91f599ede674a"
uuid = "e17b2a0c-0bdf-430a-bd0c-3a23cae4ff39"
version = "1.0.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.extxyz_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "PCRE2_jll", "Pkg", "libcleri_jll"]
git-tree-sha1 = "54c88fe9d4c57d59f3f51ac09f82e1e6ff7ef6c8"
uuid = "6ecdc6fc-93a8-5528-aee3-ac7ae1c60be7"
version = "0.1.0+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libcleri_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "PCRE2_jll", "Pkg"]
git-tree-sha1 = "9cf7bea1bc957a2133fe6ea57c95b36b5addc364"
uuid = "cdc7adba-bef8-5cba-a7ee-c792dee3081e"
version = "0.12.1+3"

[[deps.micromamba_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "f872551d05fed7616d51462992d7247f56636380"
uuid = "f8abcde7-e9b7-5caa-b8af-a437887ae8e4"
version = "1.4.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═41c2564a-ce6b-11ed-0ac3-2752d1deec22
# ╠═bf7a1e75-61f9-4255-9668-6f3a25726f1a
# ╠═7f0b1f32-3cac-4707-b521-b770b451c74f
# ╠═985dfa48-53e2-4573-bddb-fac174e2b36b
# ╟─0b718e36-3aea-4359-8414-26416277893a
# ╠═8279854f-cc71-4f95-9d22-107085ac9886
# ╠═8d151a24-3566-4eb3-af55-0d56a5701d3a
# ╠═aac6c11a-7e01-48db-8042-74fd13316ae9
# ╠═4e50d8d2-2dcb-4c0c-a4cb-93c33d89998f
# ╠═6d035aa0-48b3-4965-bbb4-7c810eedf834
# ╠═b60a0a84-b2ca-4f9e-b389-352483b1f65c
# ╟─f514edad-b5af-4a0b-99ec-d0af6b8b7423
# ╠═fea940f8-51a0-4069-9742-4e4269bf853f
# ╠═38abf557-84c9-4d6c-b7e4-a039d08a5309
# ╟─0d9e0de7-de2c-4040-ae2b-5adfbcd45987
# ╠═f7c445db-6a5c-43e6-bce0-6ea99cf11b37
# ╟─b85218e6-5d9d-48e4-a36f-88b02404bed7
# ╠═dca85534-2f5b-4ab3-b517-2853f797faaa
# ╠═5bd3263f-85a2-4fae-b3c8-62ebd5d01736
# ╠═d9c11b88-cd94-4c4d-b236-e8dad58761fe
# ╠═7b78e10f-7928-483f-ba4b-1f7151d6a360
# ╠═458c33b0-2ddc-4b2c-84a4-d5c898ad2bd8
# ╠═0f52ae1c-c683-4aaa-a1db-a570bb2b6a8f
# ╠═0977492e-c906-493b-8370-06505f2f5692
# ╠═564f902d-2d0c-4b80-81eb-03ace2cc729f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
