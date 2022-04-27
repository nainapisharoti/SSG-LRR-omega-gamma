 #include "IncompressibleTurbulenceModel.H"
 #include "transportModel.H"
 #include "addToRunTimeSelectionTable.H"
 #include "makeTurbulenceModel.H"

 #include "laminarModel.H"
 #include "RASModel.H"
 #include "LESModel.H"
 #include "ReynoldsStress.H"

 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 #define createBaseTurbulenceModel(Alpha, Rho, baseModel, BaseModel, Transport) \
                                                                                \
     namespace Foam                                                             \
     {                                                                          \
         typedef BaseModel<Transport> Transport##BaseModel;                     \
         typedef RASModel<Transport##BaseModel> RAS##Transport##BaseModel;      \
         typedef LESModel<Transport##BaseModel> LES##Transport##BaseModel;      \
     }

 createBaseTurbulenceModel
 (
     geometricOneField,
     geometricOneField,
     incompressibleTurbulenceModel,
     IncompressibleTurbulenceModel,
     transportModel
 );

 #define makeRASModel(Type)                                                     \
     makeTemplatedTurbulenceModel                                               \
     (transportModelIncompressibleTurbulenceModel, RAS, Type)

 #define makeLESModel(Type)                                                     \
     makeTemplatedTurbulenceModel                                               \
     (transportModelIncompressibleTurbulenceModel, LES, Type)

#include "ssglrr_omega.H"
makeRASModel(ssglrr_omega);
