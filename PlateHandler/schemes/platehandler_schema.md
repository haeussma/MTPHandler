```mermaid
classDiagram
    Plate *-- Well
    Plate *-- Species
    Well *-- InitCondition
    Species *-- SpeciesType
    InitCondition *-- Species
    
    class Plate {
        +integer n_rows
        +integer n_columns
        +float temperature
        +str temperature_unit
        +float ph
        +Well[0..*] wells
        +int[0..*] measured_wavelengths
        +Species[0..*] species
    }
    
    class Well {
        +float[0..*] absorption
        +float[0..*] time
        +str time_unit
        +float reaction_volume
        +string volume_unit
        +InitCondition[0..*] init_conditions
        +integer x_position
        +integer y_position
        +int wavelength
    }
    
    class Species {
        +SpeciesType type*
        +str name
    }
    
    class InitCondition {
        +Species species
        +float init_conc
        +str conc_unit
        +bool was_blanked*
    }
    
    class SpeciesType {
        << Enumeration >>
        +BUFFER
        +ENZYME
        +SUBSTRATE
        +PRODUCT
        +INHIBITOR
    }
    
```