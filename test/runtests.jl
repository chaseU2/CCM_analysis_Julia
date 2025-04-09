using CCM_analysis_Julia
using Test

using DataFrames
using DelimitedFiles

@testset "CCM Analysis with real dataset" begin
    # Pfad zu Ihrer Datendatei (anpassen falls nötig)
    data_path = "C:/Users/klbal/Desktop/Internship #2/pivoted_output_2.txt"
    
    # Überprüfen ob die Datei existiert
    if isfile(data_path)
        println("\nRunning CCM analysis with real dataset...")
        
        # Temporäres Ausgabeverzeichnis
        output_dir = mktempdir()
        
        try
            # Ihr Beispielaufruf mit allen Parametern
            final_results = run_ccm_analysis(
                data_path,
                L=100,
                E=3,
                tau=2,
                THRESHOLD=0.75,
                save_plots=true,
                save_protocol=true,
                output_dir=output_dir,
                show_plots=false  # Für automatisiertes Testen deaktiviert
            )
            
            # Ausgabe der Ergebnisse
            println("\nFinal CCM Results:")
            println(final_results)
            
            # Tests
            @test final_results isa DataFrame
            @test !isempty(final_results)  # Ergebnisse sollten nicht leer sein
            @test isdir(output_dir)
            @test any(endswith(".png", f) for f in readdir(output_dir))  # Plots vorhanden?
            @test isfile(joinpath(output_dir, "protocol.txt"))  # Protokoll vorhanden?
            
        finally
            # Aufräumen
            rm(output_dir, recursive=true)
        end
    else
        @warn "Testdatei nicht gefunden unter $data_path - Test wird übersprungen"
        @test_skip "Datei nicht verfügbar"
    end
end
