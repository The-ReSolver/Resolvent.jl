@testset "Truncation type construction      " begin
    # ensure the correct truncation type is obtained from constructor
    @test TruncatedSVD(0) isa FullSVD
    randint = rand(1:10)
    @test TruncatedSVD(randint) isa TruncatedSVD{randint}
end

@testset "Truncation performed correctly    " begin
    
end
