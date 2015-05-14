package Chemical
  block System
    parameter Integer set_n(start = 1) "Number of Compnents in Stream/System";
    constant Real pi = 3.14159265359;
    parameter Boolean use_fps = false "Give true if fps system is used";
    Real g = if use_fps then 32.174 else 9.81 "Acceleration due to gravity";
    parameter Boolean use_EnergyBal = false "Give true if Energy Balance calculations are to performed" annotation(Dialog(tab = "Energy Balance", group = "Assumptions"));
    parameter Boolean OpenSys = false "Give true if it is an open system and mass flows across the system boundary" annotation(Dialog(tab = "Energy Balance", group = "Assumptions"));
    parameter Boolean use_T = false "Give true if Temperature is used in calculations" annotation(Dialog(tab = "Energy Balance", group = "Assumptions"));
    parameter Boolean known_en = false "Give true if heat capacity of any component is known" annotation(Dialog(tab = "PPDB", group = "T Dependent"));
    parameter Boolean[set_n] known_Cv = {false} "Give true if heat capacity at constant volume is known" annotation(Dialog(tab = "PPDB", group = "T Dependent"));
    parameter Real[set_n, 4] Cv_coeff_gas = zeros(set_n, 4) "Enter coefficients of gas phase heat capacity of components at constant volume of the form A+B*T+C*T^2+D*T^3 as {A,B,C,D}" annotation(Dialog(tab = "PPDB", group = "T Dependent"));
    parameter Real[set_n, 4] Cv_coeff_liq = zeros(set_n, 4) "Enter coefficients of liquid phase heat capacity of components at constant volume of the form A+B*T+C*T^2+D*T^3 as {A,B,C,D}" annotation(Dialog(tab = "PPDB", group = "T Dependent"));
    parameter Real[set_n] Cv_coeff_aq = {0} "Enter coefficients of aqueous phase heat capacity of components at constant volume, for water give zero" annotation(Dialog(tab = "PPDB", group = "T Dependent"));
    parameter Boolean[set_n] known_Cp = {false} "Give true if heat capacity at constant pressure is known" annotation(Dialog(tab = "PPDB", group = "T Dependent"));
    parameter Real[set_n, 4] Cp_coeff_gas = zeros(set_n, 4) "Enter coefficients of gas phase heat capacity of components at constant pressure of the form A+B*T+C*T^2+D*T^3 as {A,B,C,D}" annotation(Dialog(tab = "PPDB", group = "T Dependent"));
    parameter Real[set_n, 4] Cp_coeff_liq = zeros(set_n, 4) "Enter coefficients of liquid phase heat capacity of components at constant pressure of the form A+B*T+C*T^2+D*T^3 as {A,B,C,D}" annotation(Dialog(tab = "PPDB", group = "T Dependent"));
    parameter Real[set_n] Cp_coeff_aq = {0} "Enter coefficients of aqueous phase heat capacity of components at constant pressure, for water give zero" annotation(Dialog(tab = "PPDB", group = "T Dependent"));
    parameter Real[set_n, 4] delHvap = zeros(set_n, 4) "Enter coefficients for heat of vaporization of components of the form C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2)" annotation(Dialog(tab = "PPDB", group = "T Dependent"));
    parameter Real[set_n, 3] A = zeros(set_n, 3) "Enter Antoine Coefficients [log(p*) =  A-B/(T(degC)+C)] as {{A1,B1,C1},{A2,B2,C2}}" annotation(Dialog(tab = "PPDB", group = "T Dependent"));
    parameter Real[set_n] MW "Enter Molecular weight of all the components" annotation(Dialog(tab = "PPDB", group = "Constant"));
    parameter Real[set_n, 1] pc "Enter Critical Pressure of all the components" annotation(Dialog(tab = "PPDB", group = "Constant"));
    parameter Real[set_n, 1] Tc "Enter Critical Temperature of all the components" annotation(Dialog(tab = "PPDB", group = "Constant"));
    parameter Real[set_n, 1] w "Enter Acentric Factor of all the components" annotation(Dialog(tab = "PPDB", group = "Constant"));
    parameter Real[set_n, set_n] kij "Enter Binary Interaction Parameters for all the component pairs, (3,4)th element is interaction between 3rd and 4th compound" annotation(Dialog(tab = "PPDB", group = "Constant"));
    annotation(defaultComponentName = "system", defaultComponentPrefixes = "inner", Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(lineColor = {0, 0, 255}, extent = {{-150, 150}, {150, 110}}, textString = "%name"), Line(points = {{-86, -30}, {82, -30}}), Line(points = {{-82, -68}, {-52, -30}}), Line(points = {{-48, -68}, {-18, -30}}), Line(points = {{-14, -68}, {16, -30}}), Line(points = {{22, -68}, {52, -30}}), Line(points = {{74, 84}, {74, 14}}), Polygon(fillPattern = FillPattern.Solid, points = {{60, 14}, {88, 14}, {74, -18}, {60, 14}}), Text(extent = {{16, 20}, {60, -18}}, textString = "n"), Text(extent = {{-90, 82}, {74, 50}}, textString = "defaults"), Line(points = {{-82, 14}, {-42, -20}, {2, 30}}, thickness = 0.5), Ellipse(fillColor = {255, 0, 0}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-10, 40}, {12, 18}}, endAngle = 360)}));
  end System;

  model Interfaces
    extends Modelica.Icons.InterfacesPackage;

    connector FluidPort "Interface for quasi one-dimensional fluid flow in a piping network (incompressible or compressible, one or more phases, one or more substances)"
      parameter Integer n(start = 1);
      parameter Boolean EnergyBal(start = false);
      parameter Boolean OpenSys(start = false);
      parameter Boolean use_T(start = false);
      output Real[n] n_flow;
      output Real Ek if EnergyBal "Kinetic Energy of Stream";
      output Real Ep if EnergyBal "Potential Energy of Stream";
      output Real[n] U if EnergyBal "Specific Internal Energy of Stream";
      output Real[n] H if EnergyBal and OpenSys "Enthalpy of Stream";
      output Real P(start = 101325) if EnergyBal and OpenSys "Pressure Head of Stream";
      output Real T(start = 298) if use_T;
    end FluidPort;

    connector Columnport "Fluid connector with filled, large icon to be used for vectors of FluidPorts (vector dimensions must be added after dragging)"
      extends FluidPort;
      annotation(defaultComponentName = "ports_a", Icon(coordinateSystem(extent = {{-30, -200}, {30, 200}}, preserveAspectRatio = false, initialScale = 0.2, grid = {2, 2}), graphics = {Rectangle(lineColor = {0, 127, 255}, fillColor = {255, 255, 255}, pattern = LinePattern.None, extent = {{50, -200}, {-50, 200}}), Ellipse(origin = {0, 20}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid, extent = {{-30, 180}, {30, 120}}, endAngle = 360), Ellipse(origin = {0, -20}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid, extent = {{-30, -120}, {30, -180}}, endAngle = 360)}), Diagram(coordinateSystem(extent = {{-50, -200}, {50, 200}}, preserveAspectRatio = false, initialScale = 0.2, grid = {2, 2}), graphics = {Text(origin = {0, 72}, extent = {{-75, 130}, {75, 100}}, textString = "%name"), Rectangle(lineColor = {0, 127, 255}, pattern = LinePattern.None, extent = {{25, -100}, {-25, 100}}), Ellipse(origin = {0, 78}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid, extent = {{-25, 90}, {25, 40}}, endAngle = 360), Ellipse(origin = {0, -92}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid, extent = {{-25, -40}, {25, -90}}, endAngle = 360)}));
    end Columnport;

    partial connector HeatPort "Thermal port for 1-dim. heat transfer"
      flow Real Q_flow "Heat flow rate (positive if flowing from outside into the component)";
    end HeatPort;

    connector HeatPort_In "Thermal port for 1-dim. heat transfer (filled rectangular icon)"
      extends HeatPort;
      annotation(defaultComponentName = "QIn", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {191, 0, 0}, fillColor = {191, 0, 0}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(extent = {{-50, 50}, {50, -50}}, lineColor = {191, 0, 0}, fillColor = {191, 0, 0}, fillPattern = FillPattern.Solid), Text(extent = {{-120, 120}, {100, 60}}, lineColor = {191, 0, 0}, textString = "%name")}));
    end HeatPort_In;

    connector HeatPort_Out "Thermal port for 1-dim. heat transfer (unfilled rectangular icon)"
      extends HeatPort;
      annotation(defaultComponentName = "QOut", Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(extent = {{-50, 50}, {50, -50}}, lineColor = {191, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-100, 120}, {120, 60}}, lineColor = {191, 0, 0}, textString = "%name")}), Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {191, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}));
    end HeatPort_Out;

    connector Mports_In
      extends FluidPort;
      annotation(defaultComponentName = "MportsIn", Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-50, -200}, {50, 200}}, initialScale = 0.2), graphics = {Text(extent = {{-75, 130}, {75, 100}}, textString = "%name"), Rectangle(extent = {{-25, 100}, {25, -100}}, lineColor = {0, 0, 0}), Ellipse(extent = {{-25, 90}, {25, 40}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-25, 25}, {25, -25}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-25, -40}, {25, -90}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-15, -50}, {15, -80}}, lineColor = {0, 127, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-15, 15}, {15, -15}}, lineColor = {0, 127, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-15, 50}, {15, 80}}, lineColor = {0, 127, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-50, -200}, {50, 200}}, initialScale = 0.2), graphics = {Rectangle(extent = {{-50, 200}, {50, -200}}, lineColor = {0, 127, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-50, 180}, {50, 80}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-50, 50}, {50, -50}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-50, -80}, {50, -180}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-30, 30}, {30, -30}}, lineColor = {0, 127, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-30, 100}, {30, 160}}, lineColor = {0, 127, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-30, -100}, {30, -160}}, lineColor = {0, 127, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}));
    end Mports_In;

    connector port_Out "Generic fluid connector"
      extends FluidPort;
      annotation(defaultComponentName = "portOut", Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics = {Ellipse(extent = {{-40, 40}, {40, -40}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-150, 110}, {150, 50}}, textString = "%name")}), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics = {Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 127, 255}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid)}));
    end port_Out;

    connector Mports_Out "Fluid connector with filled, large icon to be used for vectors of FluidPorts (vector dimensions must be added after dragging)"
      extends FluidPort;
      annotation(defaultComponentName = "MportsOut", Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-50, -200}, {50, 200}}, initialScale = 0.2), graphics = {Text(extent = {{-75, 130}, {75, 100}}, textString = "%name"), Rectangle(extent = {{25, -100}, {-25, 100}}, lineColor = {0, 127, 255}), Ellipse(extent = {{-25, 90}, {25, 40}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-25, 25}, {25, -25}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-25, -40}, {25, -90}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid)}), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-50, -200}, {50, 200}}, initialScale = 0.2), graphics = {Rectangle(extent = {{50, -200}, {-50, 200}}, lineColor = {0, 127, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-50, 180}, {50, 80}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-50, 50}, {50, -50}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-50, -80}, {50, -180}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid)}));
    end Mports_Out;

    connector port_In
      extends FluidPort;
      annotation(defaultComponentName = "portIn", Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics = {Ellipse(extent = {{-40, 40}, {40, -40}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-30, 30}, {30, -30}}, lineColor = {0, 127, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-150, 110}, {150, 50}}, textString = "%name")}), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics = {Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 127, 255}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-80, 80}, {80, -80}}, lineColor = {0, 127, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}));
    end port_In;
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end Interfaces;

  package Specify
    import Chemical.Utilities;
    extends Modelica.Icons.SourcesPackage;

    block point
      extends Baseclass;
      parameter Boolean known_u = false "Give true if you know the Velocity of Stream" annotation(Dialog(tab = "Energy Transfer", group = "Mechanical Energy"));
      parameter Real u = 0 "Enter Velocity of stream, if known" annotation(Dialog(tab = "Energy Transfer", group = "Mechanical Energy"));
      parameter Real z = 0 "Enter vertical postion of stream/state from reference level" annotation(Dialog(tab = "Energy Transfer", group = "Mechanical Energy"));
      parameter Boolean known_P = true "Give true if Pressure of stream is known" annotation(Dialog(tab = "State"));
      parameter Boolean known_T = true "Give true if Temperature of stream is known" annotation(Dialog(tab = "State"));
      parameter Real P = 101325 "Enter Pressure of stream/state, if known, else 1" annotation(Dialog(tab = "State"));
      parameter Real V = 0 "Enter Specific Volume of stream/state, if known, else 0" annotation(Dialog(tab = "State"));
      parameter Real T = 25 "Enter Temperature in appropriate dimension" annotation(Dialog(tab = "State"));
      parameter Real Tref = 0 "Enter Reference Temperature in appropriate dimension" annotation(Dialog(tab = "State"));
      parameter Boolean Gas = true "Give true if the stream is in gas phase, false if it is in liquid phase" annotation(Dialog(tab = "State"));
      parameter Boolean Aqueous = false "Give true if the stream is in Aqueous phase" annotation(Dialog(tab = "State"));
      parameter Boolean[system.set_n] known_U "Give true if specific Internal Energy of stream is known" annotation(Dialog(tab = "Energy Transfer", group = "Internal Energy"));
      parameter Real[system.set_n] U "Enter specific internal energy of stream/state, if known, else 0" annotation(Dialog(tab = "Energy Transfer", group = "Internal Energy"));
      parameter Boolean[system.set_n] known_H "Give true if specific Enthalpy of component(s) is/are known for eg.{true,false}" annotation(Dialog(tab = "Energy Transfer", group = "Enthalpy"));
      parameter Real[system.set_n] H "Enter specific enthalpy of stream/state, if known, else 0" annotation(Dialog(tab = "Energy Transfer", group = "Enthalpy"));
      parameter Real[system.set_n, 1] lambda = zeros(system.set_n, 1) "Enter enthalpy change(s) for component(s) due to phase change from reference state(s), else zero for eg. {{lambda1, lambda2, 0}}" annotation(Dialog(tab = "Energy Transfer", group = "Enthalpy"));
      Real v if system.use_EnergyBal;
      Real[system.set_n, 4] U1 if system.use_EnergyBal;
      Real[system.set_n, 4] H1 if system.use_EnergyBal and system.OpenSys;
      Real[system.set_n, 4] Cv if system.use_EnergyBal;
      Real[system.set_n, 4] Cp if system.use_EnergyBal and system.OpenSys;
    equation
      if system.use_EnergyBal then
        if known_n then
          if known_VF then
            sum(port.n_flow) = n_flow * rho;
          else
            sum(port.n_flow) = n_flow;
          end if;
        end if;
        if known_u then
          v = u;
        end if;
        port.Ek = 0.5 * sum(port.n_flow) * v ^ 2;
        port.Ep = sum(port.n_flow) * system.g * z;
        for i in 1:system.set_n loop
          if known_U[i] then
            port.U[i] = U[i];
            Cv[i, :] = {0, 0, 0, 0};
            U1[i, :] = {0, 0, 0, 0};
            if system.known_en then
              U1[i, :] = {0, 0, 0, 0};
            end if;
          elseif system.known_Cv[i] then
            if Gas then
              Cv[i, :] = system.Cv_coeff_gas[i, :];
            elseif Aqueous then
              Cv[i, :] = {system.Cv_coeff_aq[i], 0, 0, 0};
            else
              Cv[i, :] = system.Cv_coeff_liq[i, :];
            end if;
            U1[i, 1] = Cv[i, 1] * (port.T - Tref);
            for j in 2:4 loop
              U1[i, j] = U1[i, j - 1] + Cv[i, j] / j * (port.T ^ j - Tref ^ j);
            end for;
            port.U[i] = U1[i, 4];
          else
            U1[i, :] = {0, 0, 0, 0};
            Cv[i, :] = {0, 0, 0, 0};
          end if;
        end for;
        if system.use_T then
          if known_T then
            port.T = T;
          end if;
        end if;
        if known_P then
          port.P = P;
        end if;
        if system.OpenSys then
          for i in 1:system.set_n loop
            if known_H[i] then
              port.H[i] = H[i] + lambda[i, 1];
              H1[i, :] = {0, 0, 0, 0};
              Cp[i, :] = {0, 0, 0, 0};
            elseif system.known_Cp[i] then
              if Gas then
                Cp[i, :] = system.Cp_coeff_gas[i, :];
              elseif Aqueous then
                Cp[i, :] = {system.Cp_coeff_aq[i], 0, 0, 0};
              else
                Cp[i, :] = system.Cp_coeff_liq[i, :];
              end if;
              H1[i, 1] = Cp[i, 1] * (port.T - Tref);
              for j in 2:4 loop
                H1[i, j] = H1[i, j - 1] + Cp[i, j] / j * (port.T ^ j - Tref ^ j);
              end for;
              port.H[i] = H1[i, 4] + lambda[i, 1];
            else
              H1[i, :] = {0, 0, 0, 0};
              Cp[i, :] = {0, 0, 0, 0};
            end if;
          end for;
          for i in 1:system.set_n loop
            port.H[i] = port.U[i] + port.P * V;
          end for;
        end if;
      end if;
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Text(origin = {0, 85}, extent = {{-100, 150}, {100, 20}}, textString = "%name"), Rectangle(origin = {0, 0}, extent = {{-50, 50}, {50, -50}})}));
    end point;

    block Dry_Basis
      extends point;
      Real[system.set_n - 1] X_dry "Mole fraction on Dry basis";
      Real H2O_drystack "Mole ratio of H2O to dry stack gas";
    equation
      for i in 1:system.set_n - 1 loop
        X_dry[i] = port.n_flow[i] / sum(port.n_flow[1:system.set_n - 1]);
      end for;
      H2O_drystack = port.n_flow[system.set_n] / sum(port.n_flow[1:system.set_n - 1]);
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(origin = {0, 0.23}, fillColor = {195, 192, 177}, fillPattern = FillPattern.Solid, extent = {{-50.12, 49.88}, {50.12, -50.35}})}));
    end Dry_Basis;

    block gas_ideal
      extends point;
      parameter Real R "Enter Gas Constant in appropriate dimensions" annotation(Dialog(tab = "State", group = "Constants"));
    equation
      Port.P * V = sum(port.n_flow) * R * (port.T + 273);
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Text(origin = {0, 85}, extent = {{-100, 150}, {100, 20}}, textString = "%name"), Rectangle(origin = {7.72665, -18.0332}, fillColor = {88, 219, 221}, fillPattern = FillPattern.Solid, extent = {{-57.848, 68.1456}, {42.3913, -32.08}})}));
    end gas_ideal;

    block gas_virial
      extends point;
      parameter Real R "Enter Gas Constants in appropriate dimensions" annotation(Dialog(tab = "State", group = "Constants"));
      parameter Real Tc "Enter Critical Temperature in appropriate dimension" annotation(Dialog(tab = "State"));
      parameter Real Pc "Enter Critical Pressure in appropriate dimension" annotation(Dialog(tab = "State"));
      parameter Real PAf "Enter appropriate Pitzer Accentric Factor" annotation(Dialog(tab = "State"));
      Real Tr, Pr "Reduced Temperature and Pressure";
      Real B0, B1, B "Virial Coefficient";
      Real V1 "Molar Volume";
    equation
      Tr = T / Tc;
      Pr = Pressure / Pc;
      B0 = 0.083 - 0.422 / Tr ^ 1.6;
      B1 = 0.139 - 0.172 / Tr ^ 4.2;
      B = R * Tc / Pc * (B0 + PAf * B1);
      V1 = V / sum(port.n_flow);
      Port.P = R * T / V1 * (1 + B / V1);
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Text(origin = {0, 85}, extent = {{-100, 150}, {100, 20}}, textString = "%name"), Rectangle(origin = {7.72335, -18.03}, fillColor = {104, 170, 249}, fillPattern = FillPattern.Solid, extent = {{-57.85, 68.15000000000001}, {42.3913, -32.08}})}));
    end gas_virial;

    block gas_Raoults
      //import Chemical.Utilities;
      extends gas_ideal;
      parameter Real[3] A "Enter Antoine Coefficients [log(p*) =  A-B/(T(degC)+C)] as {A,B,C}" annotation(Dialog(tab = "State", group = "Constants"));
      Real Tdp(start = 0) "Dew point temperature of the gas";
      Real superheat(start = 0) "Degree of Superheat of the gas";
      Real pstar "Vapour Pressure of Pure component at T";
      Real[system.set_n] p "Partial Pressure of Component";
      parameter Boolean single_condensable "Give true if only one condensable species is present" annotation(Dialog(tab = "State"));
    equation
      for i in 1:system.set_n loop
        p[i] = port.n_flow[i] / sum(port.n_flow) * Port.P;
      end for;
      pstar = Utilities.Antoine_equ(A, T - 273);
      if single_condensable and not known_X[1] then
        port.n_flow[1] / sum(port.n_flow) = pstar / Port.P;
      end if;
      p[1] = Utilities.Antoine_equ(A, Tdp);
      superheat = T - 273 - Tdp;
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(origin = {-0.23, 0.23}, fillColor = {127, 191, 145}, fillPattern = FillPattern.Solid, extent = {{-49.88, 49.88}, {50.35, -50.35}})}));
    end gas_Raoults;

    block HumidAir
      extends gas_ideal;
      parameter Real[3] A "Enter Antoine Coefficients [log(p*) =  A-B/(T(degC)+C)] as {A,B,C}" annotation(Dialog(tab = "State", group = "Constants"));
      Real Tdp "Dew point of the gas";
      Real superheat "Degree of Superheat of the gas";
      Real pstar "Vapour Pressure of Pure component at T";
      parameter Real Hr "Percent Relative Humidity" annotation(Dialog(tab = "State"));
      Real[system.set_n] p "Partial Pressure of Component";
    equation
      for i in 1:system.set_n loop
        p[i] = port.n_flow[i] / sum(port.n_flow) * Port.P;
      end for;
      p[1] = Hr / 100 * pstar / 760 * 1.01325;
      pstar = Antoine_equ(A, T - 273);
      p[1] / 1.01325 * 760 = Antoine_equ(A, Tdp);
      superheat = T - 273 - Tdp;
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(origin = {0.23, -0.228384}, fillColor = {152, 189, 209}, fillPattern = FillPattern.Solid, extent = {{-50.35, 50.35}, {49.88, -49.88}})}));
    end HumidAir;

    block gas_Henrys
      extends gas_ideal;
      parameter Real[3] A "Enter Antoine Coefficients [log(p*) =  A-B/(T(degC)+C)] as {A,B,C}" annotation(Dialog(tab = "State", group = "Constants"));
      Real Tdp "Dew point of the gas";
      Real superheat "Degree of Superheat of the gas";
      Real pstar "Vapour Pressure of Pure component at T";
      Real[system.set_n] p "Partial Pressure of Component";
      parameter Boolean single_condensable "Give true if only one condensable species is present" annotation(Dialog(tab = "State"));
    equation
      for i in 1:system.set_n loop
        p[i] = port.n_flow[i] / sum(port.n_flow) * Port.P;
      end for;
      log10(pstar) = A[1] - A[2] / (T - 273 + A[3]);
      if single_condensable and not known_X[1] then
        port.n_flow[1] / sum(port.n_flow) = pstar / Port.P;
      end if;
      log10(p[1]) = A[1] - A[2] / (Tdp + A[3]);
      superheat = T - 273 - Tdp;
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(origin = {-0.236768, 0}, fillColor = {194, 149, 166}, fillPattern = FillPattern.Solid, extent = {{-49.88, 50.12}, {50.3484, -50.12}})}));
    end gas_Henrys;

    model Interm
      outer Chemical.System system;
      parameter Boolean known_n(start = false) "Give true if you know the total flow rate" annotation(Dialog(tab = "Mass Transfer"));
      parameter Boolean[system.set_n] known_X "Give true if you know the mass/mole fraction of a component, and false otherwise e.g. {true,false}" annotation(Dialog(tab = "Mass Transfer"));
      parameter Boolean[system.set_n] known_comp "Give true if you know the mass/molar flow rate of a component, and false otherwise e.g. {true,false}" annotation(Dialog(tab = "Mass Transfer"));
      parameter Real n_flow(start = 1) "Give the total mass/molar flow rate" annotation(Dialog(tab = "Mass Transfer"));
      parameter Real[system.set_n] X "Give the mass/mole fraction of the components e.g. {0.1,0.9}" annotation(Dialog(tab = "Mass Transfer"));
      parameter Real[system.set_n] comp_flow "Give the mass/molar flow rate of the components e.g. {0.1,0.9}" annotation(Dialog(tab = "Mass Transfer"));
      Chemical.Interfaces.port In(n = system.set_n) if not system.use_EnergyBal annotation(Placement(visible = true, transformation(origin = {-140, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
      Chemical.Interfaces.port Out(n = system.set_n) if not system.use_EnergyBal annotation(Placement(visible = true, transformation(origin = {120, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, -3.55271e-015}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    equation
      if known_n then
        sum(Out.n_flow) = n_flow;
      end if;
      for i in 1:system.set_n loop
        if known_X[i] then
          Out.n_flow[i] / sum(Out.n_flow) = X[i];
        elseif known_comp[i] then
          Out.n_flow[i] = comp_flow[i];
        end if;
      end for;
      Out.n_flow = In.n_flow;
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(fillColor = {81, 133, 254}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, 45}, {100, -45}}), Text(origin = {29.0398, -13.5831}, lineColor = {255, 255, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-88.76000000000001, 22.72}, {28.8068, 3.97789}}, textString = "Intermediate")}));
    end Interm;

    block pipe
      extends point;
      parameter Real ID "Enter Internal Diameter of Pipe (in m/ft), if present, else 1" annotation(Dialog(group = "Dimension"));
      Real Area;
      Real V_flow;
    equation
      Area = system.pi * (ID / 2) ^ 2;
      V_flow = v * Area;
      if not known_u then
        V_flow = sum(port.n_flow) / rho;
      end if;
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
    end pipe;

    block Baseclass
      parameter Boolean known_n = false "Give true if you know the total flow rate" annotation(Dialog(tab = "Mass Transfer"));
      parameter Boolean[system.set_n] known_X "Give true if you know the mass/mole fraction of a component, and false otherwise e.g. {true,false}" annotation(Dialog(tab = "Mass Transfer"));
      parameter Boolean[system.set_n] known_comp "Give true if you know the mass/molar flow rate of a component, and false otherwise e.g. {true,false}" annotation(Dialog(tab = "Mass Transfer"));
      parameter Real n_flow(start = 1) "Give the total mass/molar flow rate" annotation(Dialog(tab = "Mass Transfer"));
      parameter Real[system.set_n] X "Give the mass/mole fraction of the components e.g. {0.1,0.9}" annotation(Dialog(tab = "Mass Transfer"));
      parameter Real[system.set_n] comp_flow "Give the mass/molar flow rate of the components e.g. {0.1,0.9}" annotation(Dialog(tab = "Mass Transfer"));
      parameter Boolean known_VF = false "Give true if n_flow is expressed as Volumetric flow rate of stream" annotation(Dialog(tab = "Mass Transfer"));
      parameter Real rho = 1 "Enter density of fluid, if known, else 1" annotation(Dialog(tab = "Mass Transfer"));
      outer Chemical.System system;
      Chemical.Interfaces.port_Out port(n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
    equation
      if known_n and not system.use_EnergyBal then
        sum(port.n_flow) = n_flow;
      end if;
      for i in 1:system.set_n loop
        if known_X[i] then
          port.n_flow[i] / sum(port.n_flow) = X[i];
        elseif known_comp[i] then
          port.n_flow[i] = comp_flow[i];
        end if;
      end for;
    end Baseclass;

    block Soln
      extends Baseclass;
      parameter Boolean known_u = false "Give true if you know the Velocity of Stream" annotation(Dialog(tab = "Energy Transfer", group = "Mechanical Energy"));
      parameter Real u = 0 "Enter Velocity of stream, if known" annotation(Dialog(tab = "Energy Transfer", group = "Mechanical Energy"));
      parameter Real z = 0 "Enter vertical postion of stream/state from reference level" annotation(Dialog(tab = "Energy Transfer", group = "Mechanical Energy"));
      parameter Boolean known_P = true "Give true if Pressure of stream is known" annotation(Dialog(tab = "State"));
      parameter Boolean known_T = true "Give true if Temperature of stream is known" annotation(Dialog(tab = "State"));
      parameter Real P = 101325 "Enter Pressure of stream/state, if known, else 1" annotation(Dialog(tab = "State"));
      parameter Real V = 0 "Enter Specific Volume of stream/state, if known, else 0" annotation(Dialog(tab = "State"));
      parameter Real T = 25 "Enter Temperature in appropriate dimension" annotation(Dialog(tab = "State"));
      parameter Real Tref = 25 "Enter Reference Temperature in appropriate dimension" annotation(Dialog(tab = "State"));
      parameter Boolean[system.set_n] known_U "Give true if specific Internal Energy of stream is known" annotation(Dialog(tab = "Energy Transfer", group = "Internal Energy"));
      parameter Real[system.set_n] U "Enter specific internal energy of stream/state, if known, else 0" annotation(Dialog(tab = "Energy Transfer", group = "Internal Energy"));
      parameter Boolean[system.set_n] known_H "Give true if specific Enthalpy of component(s) is/are known for eg.{true,false}" annotation(Dialog(tab = "Energy Transfer", group = "Enthalpy"));
      parameter Real[system.set_n] H "Enter specific enthalpy of stream/state, if known, else 0" annotation(Dialog(tab = "Energy Transfer", group = "Enthalpy"));
      parameter Boolean[system.set_n] PhaseChange "Give true if there is phase change after an operation" annotation(Dialog(tab = "Energy Transfer", group = "Enthalpy"));
      parameter Real[system.set_n] delH "Enter specific enthalpy of solution of components" annotation(Dialog(tab = "Energy Transfer", group = "Enthalpy"));
      Real v if system.use_EnergyBal;
      Real Pressure if system.use_EnergyBal;
      Real[system.set_n, 4] U1 if system.use_EnergyBal and system.known_en;
      Real[system.set_n, 4] H1 if system.use_EnergyBal and system.OpenSys and system.known_en;
      Real[system.set_n] HSol if system.use_EnergyBal and system.OpenSys;
    equation
      if system.use_EnergyBal then
        if known_n then
          if known_VF then
            sum(port.n_flow) = n_flow * rho;
          else
            sum(port.n_flow) = n_flow;
          end if;
        end if;
        if known_u then
          v = u;
        end if;
        port.Ek = 0.5 * sum(port.n_flow) * v ^ 2;
        port.Ep = sum(port.n_flow) * system.g * z;
        for i in 1:system.set_n loop
          if known_U[i] then
            port.U[i] = U[i];
            U1[i, :] = {0, 0, 0, 0};
          elseif known_Cv[i] then
            U1[i, 1] = Cv_aq[i, 1] * (port.T - Tref);
            for j in 2:4 loop
              U1[i, j] = U1[i, j - 1] + Cv_aq[i, j] / j * (port.T ^ j - Tref ^ j);
            end for;
            port.U[i] = U1[i, 4];
          else
            U1[i, :] = {0, 0, 0, 0};
          end if;
        end for;
        if system.use_T then
          if known_T then
            port.T = T;
          end if;
        end if;
        if known_P then
          Pressure = P;
        end if;
        if system.OpenSys then
          for i in 1:system.set_n loop
            if known_H[i] then
              port.H[i] = H[i];
              H1[i, :] = {0, 0, 0, 0};
            elseif known_Cp[i] then
              H1[i, 1] = Cp_aq[i, 1] * (port.T - Tref);
              for j in 2:4 loop
                H1[i, j] = H1[i, j - 1] + Cp_aq[i, j] / j * (port.T ^ j - Tref ^ j);
              end for;
              if known_HSol[i] then
                HSol[i] = delHSol[i];
              end if;
              port.H[i] = H1[i, 4] + delHSol[i];
            else
              H1[i, :] = {0, 0, 0, 0};
            end if;
          end for;
          for i in 1:system.set_n loop
            port.H[i] = port.U[i] + Pressure * V;
          end for;
          port.P = Pressure;
        end if;
      end if;
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Text(origin = {0, 85}, extent = {{-100, 150}, {100, 20}}, textString = "%name"), Rectangle(origin = {0, 0}, extent = {{-50, 50}, {50, -50}})}));
    end Soln;
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end Specify;

  package Material_Bal
    package EOS
      //Equations of State
      extends Modelica.Icons.VariantsPackage;

      model Virial
        parameter Boolean[4] known_state "Give true if the state property is known in the oreder {P,V,n,T} for eg. {true,true,false,true}";
        parameter Real Pressure "Enter Pressure in appropriate dimension";
        parameter Real Volume "Enter Volume in appropriate dimension";
        parameter Real Moles "Enter no. of moles in appropriate dimension";
        parameter Real Temperature "Enter Temperature in appropriate dimension";
        parameter Real R "Enter Gas Constants in appropriate dimensions";
        parameter Real Tc "Enter Critical Temperature in appropriate dimension";
        parameter Real Pc "Enter Critical Pressure in appropriate dimension";
        parameter Real PAf "Enter appropriate Pitzer Accentric Factor";
        Real P, V, n, T;
        Real Tr, Pr "Reduced Temperature and Pressure";
        Real B0, B1, B "Virial Coefficient";
        Real V1 "Molar Volume";
      equation
        if known_state[1] then
          P = Pressure;
        end if;
        if known_state[2] then
          V = Volume;
        end if;
        if known_state[3] then
          n = Moles;
        end if;
        if known_state[4] then
          T = Temperature;
        end if;
        Tr = T / Tc;
        Pr = P / Pc;
        B0 = 0.083 - 0.422 / Tr ^ 1.6;
        B1 = 0.139 - 0.172 / Tr ^ 4.2;
        B = R * Tc / Pc * (B0 + PAf * B1);
        V1 = V / n;
        P = R * T / V1 * (1 + B / V1);
        annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(fillColor = {132, 179, 222}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-4, 13}, extent = {{-61.83, 27.87}, {75.4131, -47.5421}}, textString = "Virial Equation \n of State")}));
      end Virial;

      model Ideal
        parameter Boolean[4] known_state "Give true if the state property is known in the oreder {P,V,n,T} for eg. {true,true,false,true}";
        parameter Real Pressure "Enter Pressure in appropriate dimension";
        parameter Real Volume "Enter Volume in appropriate dimension";
        parameter Real Moles "Enter no. of moles in appropriate dimension";
        parameter Real Temperature "Enter Temperature in appropriate dimension";
        parameter Real R "Enter Gas Constants in appropriate dimensions";
        Real P, V, n, T;
      equation
        if known_state[1] then
          P = Pressure;
        end if;
        if known_state[2] then
          V = Volume;
        end if;
        if known_state[3] then
          n = Moles;
        end if;
        if known_state[4] then
          T = Temperature;
        end if;
        P * V = n * R * T;
        annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(fillColor = {136, 232, 229}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-4, 13}, extent = {{-61.83, 27.87}, {75.41, -47.54}}, textString = "Ideal Equation 
 of State")}));
      end Ideal;

      model SRK
        parameter Boolean[4] known_state "Give true if the state property is known in the oreder {P,V,n,T} for eg. {true,true,false,true}";
        parameter Real Pressure "Enter Pressure in appropriate dimension";
        parameter Real Volume "Enter Volume in appropriate dimension";
        parameter Real Moles "Enter no. of moles in appropriate dimension";
        parameter Real Temperature "Enter Temperature in appropriate dimension";
        parameter Real R "Enter Gas Constants in appropriate dimensions";
        parameter Real Tc "Enter Critical Temperature in appropriate dimension";
        parameter Real Pc "Enter Critical Pressure in appropriate dimension";
        parameter Real PAf "Enter appropriate Pitzer Accentric Factor";
        Real P, V, n, T;
        Real Tr "Reduced Temperature ";
        Real a, b, m, alpha;
        Real V1 "Molar Volume";
      equation
        if known_state[1] then
          P = Pressure;
        end if;
        if known_state[2] then
          V = Volume;
        end if;
        if known_state[3] then
          n = Moles;
        end if;
        if known_state[4] then
          T = Temperature;
        end if;
        Tr = T / Tc;
        a = 0.42747 * (R * Tc) ^ 2 / Pc;
        b = 0.08663999999999999 * R * Tc / Pc;
        m = 0.48508 + 1.55171 * PAf - 0.1561 * PAf ^ 2;
        alpha = (1 + m * (1 - sqrt(Tr))) ^ 2;
        V1 = V / n;
        P = R * T / (V1 - b) - alpha * a / (V1 * (V1 + b));
        annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(fillColor = {137, 235, 178}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-4, 13}, extent = {{-61.83, 27.87}, {75.41, -47.54}}, textString = "SRK Equation 
 of State")}));
      end SRK;

      model CF
        parameter Boolean[4] known_state "Give true if the state property is known in the oreder {P,V,n,T} for eg. {true,true,false,true}";
        parameter Real Pressure "Enter Pressure in appropriate dimension";
        parameter Real Volume "Enter Volume in appropriate dimension";
        parameter Real Moles "Enter no. of moles in appropriate dimension";
        parameter Real Temperature "Enter Temperature in appropriate dimension";
        parameter Real R "Enter Gas Constants in appropriate dimensions";
        parameter Real z "Enter Compressibility Factor";
        Real P, V, n, T;
        Real V1 "Molar Volume";
      equation
        if known_state[1] then
          P = Pressure;
        end if;
        if known_state[2] then
          V = Volume;
        end if;
        if known_state[3] then
          n = Moles;
        end if;
        if known_state[4] then
          T = Temperature;
        end if;
        V1 = V / n;
        z = P * V1 / (R * T);
        annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(fillColor = {100, 243, 148}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-14, 28}, extent = {{-61.83, 27.87}, {95.55, -76.11}}, textString = "Compressibility Factor \n Equation of State")}));
      end CF;
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
    end EOS;
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end Material_Bal;

  package Utilities
    extends Modelica.Icons.UtilitiesPackage;

    function Limiting_Reagent
      extends Modelica.Icons.Function;
      input Integer n;
      input Real[n] o;
      input Real[n] c;
      output Integer LR;
    protected
      Real[n] flag1;
      Real[n] flag2;
      Integer i;
      Integer nR;
    algorithm
      for i in 1:n loop
        if c[i] < 0 then
          nR := i;
        else
          break;
        end if;
      end for;
      LR := 1;
      for i in 1:nR - 1 loop
        flag1[i] := o[i + 1] / o[1];
        flag2[i] := c[i + 1] / c[1];
        if flag1[i] < flag2[i] then
          LR := i + 1;
        end if;
      end for;
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
    end Limiting_Reagent;

    function Antoine_equ
      extends Modelica.Icons.Function;
      input Real[3] A;
      input Real T;
      output Real pstar;
    algorithm
      pstar := 10 ^ (A[1] - A[2] / (T + A[3]));
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
    end Antoine_equ;

    function VLE
      //T = [n x y P T_guess dT]
      extends Modelica.Icons.Function;
      input Integer n;
      input Real[n] x;
      input Real[n] y;
      input Real P;
      input Real T_guess;
      input Real dT;
      input Real[n, 3] A;
      output Real T;
      output Real[n] Ps;
    protected
      Integer flag1 = -1;
      Integer flag2 = -1;
    algorithm
      while 1 > 0 loop
        for i in 1:n loop
          Ps[i] := Antoine_equ(A[i, :], T_guess);
          y[i] := x[i] * Ps[i] / P;
        end for;
        if abs(sum(y) - 1) <= 0.0001 then
          T := T_guess;
          break;
        end if;
        if sum(y) < 1 then
          if flag1 > 0 then
            dT := dT / 2;
          end if;
          T_guess := T_guess + dT;
          flag2 := 1;
        elseif sum(y) > 1 then
          if flag2 > 0 then
            dT := dT / 2;
          end if;
          T_guess := T_guess - dT;
          flag1 := 1;
        end if;
      end while;
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
    end VLE;

    function PR_EOS
      extends Modelica.Icons.Function;
      // Peng-Robinson EOS computation
      // Inputs:
      // z   = number of moles of all species (n x 1 column vector) [mol]
      // p   = pressure [Pa]
      // T   = temperature [K]
      // pc  = critical pressure of all components (n x 1 column vector)[Pa]
      // Tc  = critical temperature of all components (n x 1 column vector)[K]
      // w   = acentric factor of all components (n x 1 column vector)
      // k   = binary parameters (n x n symmetric matrix)
      // cp  = ideal gas heat capacity coefficients (n x m symmetric matrix) [J/mol-K]
      // DHf = standard enthalpy of formation for ideal gas (298.15 K and 1 atm) (n x 1 column vector) [J/mol]
      // ph  = phase (L or V)
      // Outputs:
      // V   = molar volume [m3/mol]
      // Z   = compressibility factor
      // phi = fugacity coefficient
      // H   = enthalpy [J/mol]
      //phi = PR_EOS(n,z,p,T,pc,Tc,w,k,ph)
      input Integer n;
      input Real[n, 1] z;
      input Real p;
      input Real T;
      input Real[n, 1] pc;
      input Real[n, 1] Tc;
      input Real[n, 1] w;
      input Real[n, n] k;
      //input Real[n,5] cpig;
      //input Real[n,1] DHf;
      input String ph;
      //output Real V;
      //output Real Z;
      output Real[n, 1] phi;
      //output Real H;
    protected
      Real V;
      Real Z;
      //Real H;
      Real R;
      constant Real e = 1 - sqrt(2);
      constant Real s = 1 + sqrt(2);
      //Tr,Hr;
      Real[n, 1] m;
      Real[n, 1] alfa;
      Real ai[n, 1];
      Real bi[n, 1];
      Real Q[n, n];
      Real a[1, 1];
      Real b[1, 1];
      Real[n, n] dQdT;
      Real[1, 1] dadT;
      Real c[4];
      Real r[3, 2];
      Real[n, 1] abar;
      Real[n, 1] bbar;
      //Real Hig[n];
      Real[n, 1] J;
      Real[n, n] J1;
      Real[n, 1] J2;
      Real[n, 1] J3;
      Real[n, 1] J4;
      Real[n, 4] J6;
      Real[3] r1;
      Real[1, n] J7;
      Real[1, n] J8;
      Real[n, 1] J9;
    algorithm
      R := 8.314;
      // m3 Pa/(mol K) = J/mol-K
      //Tr:=298.15;
      // K
      m := 0.37464 * ones(n, 1) + 1.54226 * w - 0.26992 * w .^ 2;
      alfa := (ones(n, 1) + m .* (ones(n, 1) - (T ./ Tc) .^ 0.5)) .^ 2;
      J := Tc .^ 2;
      ai := 0.45724 * R ^ 2 * J ./ pc .* alfa;
      bi := 0.07779999999999999 .* R .* Tc ./ pc;
      Q := (ai * transpose(ai)) .^ 0.5 .* (ones(n, n) - k);
      a := transpose(z) * Q * z;
      b := transpose(z) * bi;
      J1 := (pc * transpose(pc)) .^ 0.5;
      J2 := alfa .^ 0.5;
      J3 := Tc .^ 0.5;
      J4 := m ./ J3;
      J6 := Tc * transpose(Tc);
      dQdT := 0.45724 * R ^ 2 * (k - ones(n, n)) .* J6 ./ J1 .* 1 / (2 * T ^ 0.5) .* (J4 * transpose(J2) + J2 * transpose(J4));
      dadT := transpose(z) * (Q - T * dQdT) * z;
      // Coefficients of the EoS model equation
      c[1] := 1;
      c[2] := b[1, 1] - R * T / p;
      c[3] := (-3 * b[1, 1] ^ 2) - 2 * R * T / p * b[1, 1] + a[1, 1] / p;
      c[4] := b[1, 1] ^ 3 + R * T / p * b[1, 1] ^ 2 - a[1, 1] * b[1, 1] / p;
      // Roots
      r := Modelica.Math.Vectors.Utilities.roots(c);
      //r := {{0.0052, 0}, {0.0007, 0}, {0.00013715, 0}};
      for i in 1:3 loop
        r1[i] := r[i, 1];
      end for;
      if ph == "L" then
        V := min(r1);
      else
        V := max(r1);
      end if;
      Z := p * V / (R * T);
      J7 := transpose(z) * Q;
      J8 := a[1, 1] .* ones(1, n);
      J9 := 2 .* J7 - J8;
      abar := transpose(J9);
      bbar := bi;
      for i in 1:n loop
        phi[i, 1] := exp((Z - 1) * bbar[i, 1] / b[1, 1] - log((V - b[1, 1]) * Z / V) + a[1, 1] / (b[1, 1] * R * T) / (e - s) * log((V + s * b[1, 1]) / (V + e * b[1, 1])) * (1 + abar[i, 1] / a[1, 1] - bbar[i, 1] / b[1, 1]));
      end for;
      //Hig[i]:=z[i,1] * (DHf[i,1] + cpig[i,1] * (T - Tr) + cpig[i,2] .* cpig[i,3] .* (1 / tanh(cpig[i,3] / T) - 1 / tanh(cpig[i,3] / Tr)) + cpig[i,4] .* cpig[i,5] .* (1 / tanh(cpig[i,5] / T) - 1 / tanh(cpig[i,5] / Tr)));
      //HR:=R * T * (Z - 1) - dadT[1,1] / (2 * (s - 1) * b[1,1]) * log((V + s * b[1,1]) / (V + e * b[1,1]));
      //H:=sum(Hig) + HR;
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
    end PR_EOS;

    function PTFlash
      extends Modelica.Icons.Function;
      //(T,K) = [n x y P T_guess dT pc Tc w k]
      input Integer n;
      input Real[n, 1] x;
      input Real[n, 1] y;
      input Real P;
      input Real T_guess;
      input Real dT;
      input Real[n, 1] pc;
      input Real[n, 1] Tc;
      input Real[n, 1] w;
      input Real[n, n] k;
      input String ph;
      output Real T;
      output Real[n, 1] K;
    protected
      Real dT1;
      Real T_guess1;
      Real[n, 1] phiL;
      Real[n, 1] phiV;
      Integer flag1 = -1;
      Integer flag2 = -1;
    algorithm
      T_guess1 := T_guess;
      dT1 := dT;
      while 1 > 0 loop
        phiL := PR_EOS(n, x, P, T_guess1, pc, Tc, w, k, "L");
        phiV := PR_EOS(n, y, P, T_guess1, pc, Tc, w, k, "V");
        K := phiL ./ phiV;
        y := K .* x;
        if abs(sum(y) - 1) <= 0.0001 then
          T := T_guess1;
          break;
        end if;
        if sum(y) < 1 then
          if flag1 > 0 then
            dT1 := dT1 / 2;
          end if;
          T_guess1 := T_guess1 + dT1;
          flag2 := 1;
        elseif sum(y) > 1 then
          if flag2 > 0 then
            dT1 := dT1 / 2;
          end if;
          T_guess1 := T_guess1 - dT1;
          flag1 := 1;
        end if;
      end while;
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
    end PTFlash;

    function VLE_PR
      extends Modelica.Icons.Function;
      //(T,K) = [n z P T_guess dT pc Tc w k]
      input Integer n;
      input Real[n, 1] z;
      input Real P;
      input Real T_guess;
      input Real dT;
      input Real[n, 1] pc;
      input Real[n, 1] Tc;
      input Real[n, 1] w;
      input Real[n, n] k;
      output Real T;
      output Real[n, 1] K;
    protected
      Real dT1;
      Real T_guess1;
      Real[n, 1] x;
      Real[n, 1] y;
      Real[n, 1] phiL;
      Real[n, 1] phiV;
      Integer flag1 = -1;
      Integer flag2 = -1;
    algorithm
      T_guess1 := T_guess;
      dT1 := dT;
      while 1 > 0 loop
        x := z - y;
        phiL := PR_EOS(n, x, P, T_guess1, pc, Tc, w, k, "L");
        phiV := PR_EOS(n, y, P, T_guess1, pc, Tc, w, k, "V");
        K := phiL ./ phiV;
        y := K .* x;
        if abs(sum(y) - 1) <= 0.0001 then
          T := T_guess1;
          break;
        end if;
        if sum(y) < 1 then
          if flag1 > 0 then
            dT1 := dT1 / 2;
          end if;
          T_guess1 := T_guess1 + dT1;
          flag2 := 1;
        elseif sum(y) > 1 then
          if flag2 > 0 then
            dT1 := dT1 / 2;
          end if;
          T_guess1 := T_guess1 - dT1;
          flag1 := 1;
        end if;
      end while;
    end VLE_PR;
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end Utilities;

  model GenericModel
    //  parameter Integer nin "Number of inputs" annotation(Dialog(group = "Number of Ports"));
    //  parameter Integer nout "Number of outputs" annotation(Dialog(group = "Number of Ports"));
    Real[nin, system.set_n] l(start = ones(nin, system.set_n)) "for calculation only. Please ignore";
    Real[nout, system.set_n] m(start = ones(nout, system.set_n)) "for calculation only. Please ignore";
    Real[system.set_n] o "for calculation only. Please ignore";
    Real[system.set_n] p "for calculation only. Please ignore";
    Real[nin] q "for calculation only. Please ignore";
    Real[nout] r "for calculation only. Please ignore";
    Real delEk if system.use_EnergyBal "Change in Kinetic Energy";
    Real delEp if system.use_EnergyBal "Change in Potential Energy";
    Real delU if system.use_EnergyBal "Change in Internal Energy of system";
    Real delH if system.use_EnergyBal and system.OpenSys and not Bernoulli "Change in Enthalpy of system";
    Real delVHead if system.use_EnergyBal and Bernoulli "Velocity Head";
    Real delSHead if system.use_EnergyBal and Bernoulli "Static Head";
    Real delPHead if system.use_EnergyBal and system.OpenSys "Change in Pressure Head of system";
    Real[nin] Ekin if system.use_EnergyBal;
    Real[nin] Epin if system.use_EnergyBal;
    Real[nin] Uin if system.use_EnergyBal;
    Real[nin] PHeadin if system.use_EnergyBal and system.OpenSys;
    Real[nout] Ekout if system.use_EnergyBal;
    Real[nout] Epout if system.use_EnergyBal;
    Real[nout] Uout if system.use_EnergyBal;
    Real[nout] PHeadout if system.use_EnergyBal and system.OpenSys;
    Real[nin, system.set_n] H1 if system.use_EnergyBal and system.OpenSys and not Bernoulli "for calculation only. Please ignore";
    Real[nout, system.set_n] H2 if system.use_EnergyBal and system.OpenSys and not Bernoulli "for calculation only. Please ignore";
    Real[nin, system.set_n] U1 if system.use_EnergyBal "for calculation only. Please ignore";
    Real[nout, system.set_n] U2 if system.use_EnergyBal "for calculation only. Please ignore";
    Real[nin] Hin if system.use_EnergyBal and system.OpenSys and not Bernoulli;
    Real[nout] Hout if system.use_EnergyBal and system.OpenSys and not Bernoulli;
    // Real[system.set_n] Hp if system.use_EnergyBal and system.OpenSys and PhaseChange;
    outer Chemical.System system;
  equation
    for i in 1:nin loop
      for j in 1:system.set_n loop
        l[i, j] = In[i].n_flow[j];
      end for;
      q[i] = sum(l[i, :]);
    end for;
    for i in 1:nout loop
      for j in 1:system.set_n loop
        m[i, j] = Out[i].n_flow[j];
      end for;
      r[i] = sum(m[i, :]);
    end for;
    for k in 1:system.set_n loop
      o[k] = sum(l[:, k]);
      p[k] = sum(m[:, k]);
    end for;
    if system.use_EnergyBal then
      for i in 1:nin loop
        Ekin[i] = In[i].Ek;
        Epin[i] = In[i].Ep;
      end for;
      for i in 1:nout loop
        Ekout[i] = Out[i].Ek;
        Epout[i] = Out[i].Ep;
      end for;
      for i in 1:nin loop
        for j in 1:system.set_n loop
          U1[i, j] = In[i].n_flow[j] * In[i].U[j];
        end for;
        Uin[i] = sum(U1[i, :]);
      end for;
      for i in 1:nout loop
        for j in 1:system.set_n loop
          U2[i, j] = Out[i].n_flow[j] * Out[i].U[j];
        end for;
        Uout[i] = sum(U2[i, :]);
      end for;
      delEk = sum(Ekout) - sum(Ekin);
      delEp = sum(Epout) - sum(Epin);
      delU = sum(Uout) - sum(Uin);
      if system.OpenSys then
        for i in 1:nin loop
          PHeadin[i] = In[i].P;
        end for;
        for i in 1:nout loop
          PHeadout[i] = Out[i].P;
        end for;
        delPHead = sum(PHeadout) - sum(PHeadin);
        if Bernoulli then
          delVHead = sum(Ekout ./ r) - sum(Ekin ./ q);
          delSHead = sum(Epout ./ r) - sum(Epin ./ q);
        else
          for i in 1:nin loop
            for j in 1:system.set_n loop
              H1[i, j] = In[i].n_flow[j] * In[i].H[j];
            end for;
            Hin[i] = sum(H1[i, :]);
          end for;
          for i in 1:nout loop
            for j in 1:system.set_n loop
              //if PhaseChange then
              //if known_delHp[j] then
              //Hp[j] = delHp[j];
              //end if;
              // H2[i,j] = Out[i].n_flow[j] * (Out[i].H[j] + Hp[j]);
              //else
              H2[i, j] = Out[i].n_flow[j] * Out[i].H[j];
              //end if;
            end for;
            Hout[i] = sum(H2[i, :]);
          end for;
          delH = sum(Hout) - sum(Hin);
        end if;
      end if;
    end if;
    annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end GenericModel;

  package NonReactive_Sys
    extends Modelica.Icons.VariantsPackage;

    model GeneralBal
      extends Chemical.GenericModel;
      parameter Boolean known_Q "Give true if the heat transferred to/from the system is known" annotation(Dialog(tab = "Energy Transfer", group = "Heat"));
      parameter Real Heatflow "Give Heat transf. to the system as +ve else if the system is adiabatic = 0" annotation(Dialog(tab = "Energy Transfer", group = "Heat"));
      parameter Boolean Temp_change "Give true if there is temperature change" annotation(Dialog(tab = "Energy Transfer", group = "Heat"));
      parameter Boolean known_W "Give true if the work done on/by the system is known" annotation(Dialog(tab = "Energy Transfer", group = "Work"));
      parameter Real Work "Give Work done by the system as +ve else if there are no moving parts = 0" annotation(Dialog(tab = "Energy Transfer", group = "Work"));
      parameter Boolean Bernoulli = false "Give true if Heat changes are negligible" annotation(Dialog(tab = "Energy Transfer", group = "Mechanical"));
      parameter Boolean known_FLoss = false "Give true if the friction loss from the system is known" annotation(Dialog(tab = "Energy Transfer", group = "Mechanical"));
      parameter Real F_loss = 0 "Give Friction Losses from the system" annotation(Dialog(tab = "Energy Transfer", group = "Mechanical"));
      Real Q if system.use_EnergyBal and not Bernoulli "Heat transferred to/from the system from/to the surroundings";
      Real W if system.use_EnergyBal "Work done on/by the system";
      Real F if system.use_EnergyBal and Bernoulli "Friction losses from the system";
    equation
      // Mass Balance
      o = p;
      // Energy Balance
      if system.use_EnergyBal then
        if known_W then
          W = Work;
        end if;
        if Bernoulli then
          if known_FLoss then
            F = F_loss;
          end if;
          delPHead + delVHead + delSHead + F = -W / sum(p);
        else
          if not Temp_change then
            delU = 0;
          end if;
          if known_Q then
            Q = Heatflow;
          end if;
          if system.OpenSys then
            delH + delEk + delEp = Q - W;
          else
            delU + delEk + delEp = Q - W;
          end if;
        end if;
      end if;
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
    end GeneralBal;

    model Liquid_Vapour_System
      parameter Integer nin(start = 1) "Number of inputs" annotation(Dialog(group = "Number of Ports"));
      outer Chemical.System system;
      parameter Boolean known_Liquid "Give true if the composition of Liquid is known";
      parameter Boolean known_P "Give true if Pressure of the system is known is known" annotation(Dialog(tab = "State"));
      parameter Real[system.set_n - NC, 3] A "Enter Antoine Coefficients [log(p*) =  A-B/(T(degC)+C)] as {A,B,C}";
      parameter Real P "Give Total Pressure" annotation(Dialog(tab = "State"));
      parameter Real T "Give Temperature of System if known" annotation(Dialog(tab = "State"));
      parameter Integer NC(start = system.set_n) "Give number of Non-Condensable elements";
      Real Tbp "Bubble Point temperature";
      Real Tdp "Dew Point temperature";
      Real T1 "Only for calculation";
      Real P1 "Only for calculation";
      Real[system.set_n - NC] pstar "Vapour Pressure of Pure components";
      Real[system.set_n - NC] p "Partial Pressure of components";
      Chemical.Interfaces.Mports ports_a annotation(Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-5, -20}, {5, 20}}, rotation = 0), iconTransformation(origin = {-99.5, -1}, extent = {{-12.5, -50}, {12.5, 50}}, rotation = 0)));
    equation
      for i in 1:system.set_n - NC loop
        pstar[i] = Antoine_equ(A[i, :], T1 - 273);
      end for;
      if known_Liquid then
        Tbp = T1;
        for i in 1:system.set_n - NC loop
          p[i] = u[2].n_flow[i] * pstar[i];
          u[1].n_flow[i] = p[i] / P;
        end for;
        Tdp = 0;
        P1 = P;
      else
        if not known_P then
          T = T1;
          Tdp = 0;
        else
          Tdp = T1;
          P1 = P;
        end if;
        for i in 1:system.set_n - NC loop
          p[i] = u[1].n_flow[i] * P1;
          u[2].n_flow[i] = p[i] / pstar[i];
        end for;
        Tbp = 0;
      end if;
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(fillColor = {18, 69, 253}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-9.65, 10}, extent = {{-50, 50}, {77.56180000000001, -66.6078}}, textString = "Liquid-Vapour \n System")}));
    end Liquid_Vapour_System;

    model ColligativePr
      parameter Real Xs "Give Mass of Solute present in 100gm of solvent";
      parameter Real T "Give Temperature of mixture";
      parameter Boolean known_BP "Give true if Boiling point of the solution is known";
      parameter Real BP_obs "Give observed Boiling Point of Solution";
      parameter Real delHv "Give Heat of Vaporiztion of pure solvent";
      parameter Boolean known_MP "Give true if Melting point of the solution is known";
      parameter Real MP_obs "Give observed Melting Point of Solution";
      parameter Real delHm "Give Heat of fusion of pure solvent";
      parameter Real[3] A "Enter Antoine Coefficients [log(p*) =  A-B/(T(degC)+C)] as {A,B,C}";
      Real MWt "Molecular Weight of Solute";
      Real x "Solute Mole Fraction";
      Real delTb "Elevation in Boiling Point";
      Real delTm "Depression in Freezing point";
      Real pstar "Pure Solvent Vapour Pressure";
      Real pse "Effetive Solvent Vapour Pressure";
      Real BP;
      Real MP;
      constant Real R = 8.314 "Gas Constant";
    equation
      if known_BP then
        BP = BP_obs;
      elseif known_MP then
        MP = MP_obs;
      end if;
      delTb = R * 373.16 ^ 2 * x / delHv;
      delTb = BP - 373.16;
      x = Xs / MWt / (Xs / MWt + 100 / 18.06);
      pstar = Antoine_equ(A, T - 273);
      pse = (1 - x) * pstar;
      delTm = R * 273.16 ^ 2 * x / delHm;
      delTm = 273.16 - MP;
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(fillColor = {19, 212, 90}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-53.3, 27.59}, {53.3, -27.59}}, textString = "Colligative Properties 
 of Solutions")}));
    end ColligativePr;

    model LLE
      extends GeneralBal;
      parameter Real K "Give Distribution Coefficient for the system";
      Real Per_tr "Percentage transferred into Solvent";
    equation
      y[2].n_flow[1] / sum(y[2].n_flow) / (y[1].n_flow[1] / sum(y[1].n_flow)) = K;
      Per_tr = y[2].n_flow[1] / u[1].n_flow[1];
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(fillColor = {131, 64, 245}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-53.07, 29.95}, {53.07, -29.95}}, textString = "LLE")}));
    end LLE;

    model Solution
      extends Chemical.GenericModel;
      parameter Boolean Temp_change "Give true if there is a Temperature Change";
      parameter Real xnew "Give solubility at new temperature";
      parameter Real MWt_Salt "Give Molecular Weight of Anhydrous Salt";
      parameter Boolean Salt_Hydra "Give true if the Crystals are hydrated";
      parameter Real Water_Hydra "Give Water of Hydration/no. of water molecules attached with crystal";
      Real MWt_Cr "Molecular Weight of hydrated Crystal";
      Real dummy "Only for Calculation";
    equation
      if Salt_Hydra then
        MWt_Cr = MWt_Salt + Water_Hydra * 18;
        dummy = MWt_Salt / MWt_Cr;
        o[1] + o[3] = p[1] + p[3] * dummy;
        o[2] = p[2] + p[3] * (1 - dummy);
      else
        o[2] = p[2];
        o[1] + o[3] = p[1] + p[3];
        MWt_Cr = 0;
        dummy = 0;
      end if;
      if Temp_change then
        p[1] / (p[1] + p[2]) = xnew;
      end if;
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(fillColor = {70, 152, 239}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-53.3, 27.59}, {53.3, -27.59}}, textString = "Balance \n on \n Solutions")}));
    end Solution;
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end NonReactive_Sys;

  package UnitOp
    model Mixer
      extends Chemical.Icons.Mixer;
      parameter Integer nin "Number of input(s)" annotation(Dialog(group = "Number of Ports"));
      parameter Integer nout = 1 "Number of Output" annotation(Dialog(group = "Number of Ports", enable = false));
      extends Chemical.NonReactive_Sys.GeneralBal(final known_W = true, final Work = 0, final Bernoulli = false, final known_FLoss = false, final F_loss = 0);
      parameter Integer P_downstream = 1 "1:Average, 2: Maximum, 3: Minimum " annotation(Dialog(tab = "Specifications"));
      Real[nin] dummy_sum;
      Interfaces.Mports_In In[nin](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-5, -20}, {5, 20}}, rotation = 0), iconTransformation(origin = {-99.5, 1}, extent = {{-13.5, -54}, {13.5, 54}}, rotation = 0)));
      Interfaces.port_Out Out[nout](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {80, 3.33067e-015}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    equation
      for i in 1:nin loop
        dummy_sum[i] = In[i].P;
      end for;
      if P_downstream == 1 then
        Out[1].P = sum(dummy_sum) / nin;
      elseif P_downstream == 2 then
        Out[1].P = max(dummy_sum);
      elseif P_downstream == 3 then
        Out[1].P = min(dummy_sum);
      end if;
    end Mixer;

    model Flash2
      extends Icons.Column;
      import Chemical.Utilities.*;
      parameter Integer nin = 1 "Number of inputs" annotation(Dialog(group = "Number of Ports", enable = false));
      parameter Integer nout = 2 "Number of Output" annotation(Dialog(group = "Number of Ports", enable = false));
      extends Chemical.NonReactive_Sys.GeneralBal(final known_W = true, final Work = 0, final Bernoulli = false, final known_FLoss = false, final F_loss = 0);
      Chemical.Interfaces.Mports_In In[nin](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {-99.5, -2}, extent = {{-12.5, -50}, {12.5, 50}}, rotation = 0), iconTransformation(origin = {-50, 0}, extent = {{-5, -20}, {5, 20}}, rotation = 0)));
      parameter Boolean known_P "Give true if operating Pressure is known" annotation(Dialog(tab = "Specifications"));
      parameter Boolean known_T "Give true if operating Temperature is known" annotation(Dialog(tab = "Specifications"));
      parameter Real Pressure "Enter operating Pressure" annotation(Dialog(tab = "Specifications"));
      parameter Real Temperature "Enter operating Temperature" annotation(Dialog(tab = "Specifications"));
      parameter Real T_guess = 25 "Enter guess valur for equilibrium Temperature if operating temperature is not known" annotation(Dialog(tab = "Specifications"));
      Real P;
      Real T if system.use_T;
      Real[nin] dummy_sum;
      constant Real dT = 20;
      Real[system.set_n, 1] y;
      Real[system.set_n, 1] x;
      Real[system.set_n, 1] phiL if known_T;
      Real[system.set_n, 1] phiV if known_T;
      Real[system.set_n, 1] K;
      Chemical.Interfaces.Columnport Out[nout](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {0, 0}, extent = {{-23.2812, -93.125}, {23.2812, 93.125}}, rotation = 0), iconTransformation(origin = {0.5, 3.55271e-015}, extent = {{-13.5, -90}, {13.5, 90}}, rotation = 0)));
    equation
      x[:, 1] = Out[2].n_flow ./ sum(Out[2].n_flow);
      y[:, 1] = Out[1].n_flow ./ sum(Out[1].n_flow);
      if known_P then
        P = Pressure;
      else
        for i in 1:nin loop
          dummy_sum[i] = In[i].P;
        end for;
        P = sum(dummy_sum) / nin;
      end if;
      Out[1].P = P;
      Out[2].P = P;
      if system.use_T then
        if known_T then
          T = Temperature;
          phiL = PR_EOS(system.set_n, x, P, T, system.pc, system.Tc, system.w, system.kij, "L");
          phiV = PR_EOS(system.set_n, y, P, T, system.pc, system.Tc, system.w, system.kij, "V");
          K = phiL ./ phiV;
          y = x .* K;
        else
          (T, K) = PTFlash(system.set_n, x, P, T_guess, dT, system.pc, system.Tc, system.w, system.kij);
          y = x .* K;
        end if;
        for i in 1:nout - 1 loop
          Out[i].T = T;
        end for;
      end if;
      annotation(Icon(coordinateSystem(extent = {{-55, -90}, {50, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
    end Flash2;

    model FSplit
      extends Chemical.Icons.Splitter;
      extends Chemical.NonReactive_Sys.GeneralBal(final known_W = true, final Work = 0, final Bernoulli = false, final known_FLoss = false, final F_loss = 0);
      parameter Integer nin "Number of input(s)" annotation(Dialog(group = "Number of Ports"));
      parameter Integer nout "Number of output(s)" annotation(Dialog(group = "Number of Ports"));
      parameter Boolean[nout] known_Fsplit "Give true if fraction of combined inlet flow is/are known for outlet stream(s), give false for atleast one; for eg. for 3 outlets: {true,true,false}" annotation(Dialog(tab = "Product"));
      parameter Real[nout] alpha "Enter the split fractions for the outlet streams" annotation(Dialog(tab = "Product"));
      Real P;
      Real[nin] dummy_sum;
      Chemical.Interfaces.Mports_Out Out[nout](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {100, 0}, extent = {{-12.5, -50}, {12.5, 50}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -40}, {10, 40}}, rotation = 0)));
      Chemical.Interfaces.Mports_In In[nin](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {-100, 20}, extent = {{-5, -20}, {5, 20}}, rotation = 0), iconTransformation(origin = {-95, 0}, extent = {{-5, -20}, {5, 20}}, rotation = 0)));
    equation
      for i in 1:nin loop
        dummy_sum[i] = In[i].P;
      end for;
      P = sum(dummy_sum) / nin;
      for i in 1:nout loop
        Out[i].P = P;
        if known_Fsplit[i] then
          sum(Out[i].n_flow) = alpha[i] * sum(o);
        end if;
        if i > 1 then
          Out[i].T = Out[1].T;
          for j in 1:system.set_n - 1 loop
            Out[i].n_flow[j] / sum(Out[i].n_flow) = Out[1].n_flow[j] / sum(Out[1].n_flow);
          end for;
        end if;
      end for;
    end FSplit;

    model Sep
      extends Icons.Column;
      import Chemical.Utilities.*;
      extends Chemical.NonReactive_Sys.GeneralBal(final known_W = true, final Work = 0, final Bernoulli = false, final known_FLoss = false, final F_loss = 0);
      parameter Integer nin = 1 "Number of inputs" annotation(Dialog(group = "Number of Ports"));
      parameter Integer nout = 2 "Number of Output" annotation(Dialog(group = "Number of Ports"));
      parameter Boolean[nout, system.set_n] known_Ssplit "Give true if component split fraction(s) is/are known for the outlet stream(s), give false for atleast one outlet for a particular component, for eg. 2 outlets and 2 components {{true,false},{false,true}}" annotation(Dialog(tab = "Product Flash"));
      parameter Real[nout, system.set_n] beta "Enter the split fractions for the components od outlet streams" annotation(Dialog(tab = "Product"));
      Chemical.Interfaces.Mports_In In[nin](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {-99.5, -2}, extent = {{-12.5, -50}, {12.5, 50}}, rotation = 0), iconTransformation(origin = {-50, 0}, extent = {{-5, -20}, {5, 20}}, rotation = 0)));
      Chemical.Interfaces.Mports_Out Out[nout](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {100.5, -2}, extent = {{-12.5, -50}, {12.5, 50}}, rotation = 0), iconTransformation(origin = {50, 0}, extent = {{-5, -20}, {5, 20}}, rotation = 0)));
    equation
      for i in 1:nout loop
        for j in 1:system.set_n loop
          if known_Ssplit[i, j] then
            Out[i].n_flow[j] = beta[i, j] * o[j];
          end if;
        end for;
      end for;
      if known_P then
        P = Pressure;
      end if;
      annotation(Icon(coordinateSystem(extent = {{-55, -90}, {55, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Line(points = {{-50, 60}, {50, -60}}, smooth = Smooth.Bezier), Line(points = {{-50, -60}, {50, 60}}, smooth = Smooth.Bezier)}));
    end Sep;

    model Sep2
      extends Icons.Column;
      import Chemical.Utilities.*;
      extends Chemical.NonReactive_Sys.GeneralBal(final known_W = true, final Work = 0, final Bernoulli = false, final known_FLoss = false, final F_loss = 0);
      Chemical.Interfaces.Mports_In In[nin](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {-99.5, -2}, extent = {{-12.5, -50}, {12.5, 50}}, rotation = 0), iconTransformation(origin = {-50, 0}, extent = {{-5, -20}, {5, 20}}, rotation = 0)));
      parameter Integer nin = 1 "Number of inputs" annotation(Dialog(group = "Number of Ports"));
      parameter Integer nout = 2 "Number of Output" annotation(Dialog(group = "Number of Ports", enable = false));
      parameter Boolean[nout] known_Fsplit "Give true if fraction of combined inlet flow is/are known for outlet stream(s), give false for atleast one; for eg. for 3 outlets: {true,true,false}" annotation(Dialog(tab = "Product"));
      parameter Real[nout] alpha "Enter the split fractions for the outlet streams" annotation(Dialog(tab = "Product"));
      parameter Boolean[nout, system.set_n] known_Ssplit "Give true if component split fraction(s) is/are known for the outlet stream(s), give false for atleast one outlet for a particular component, for eg. 2 outlets and 2 components {{true,false},{false,true}}" annotation(Dialog(tab = "Product"));
      parameter Real[nout, system.set_n] beta "Enter the split fractions for the components of outlet streams" annotation(Dialog(tab = "Product"));
      parameter Boolean known_P "Give true if operating Pressure is known" annotation(Dialog(tab = "Specifications"));
      parameter Boolean known_T "Give true if operating Temperature is known" annotation(Dialog(tab = "Specifications"));
      parameter Real Pressure "Enter operating Pressure" annotation(Dialog(tab = "Specifications"));
      parameter Real Temperature "Enter operating Temperature" annotation(Dialog(tab = "Specifications"));
      parameter Real T_guess "Enter guess valur for equilibrium Temperature" annotation(Dialog(tab = "Specifications"));
      Real P;
      Real T if system.use_T;
      Real[nin] dummy_sum;
      constant Real dT = 20;
      Real[system.set_n, 1] y;
      Real[system.set_n, 1] x;
      Real[system.set_n, 1] phiL;
      Real[system.set_n, 1] phiV;
      Real[system.set_n, 1] K;
      Chemical.Interfaces.Columnport Out[nout](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {0, 0}, extent = {{-23.2812, -93.125}, {23.2812, 93.125}}, rotation = 0), iconTransformation(origin = {0.5, 3.55271e-015}, extent = {{-13.5, -90}, {13.5, 90}}, rotation = 0)));
    equation
      for i in 1:nout loop
        if known_Fsplit[i] then
          sum(Out[i].n_flow) = alpha[i] * sum(o);
        end if;
        for j in 1:system.set_n loop
          if known_Ssplit[i, j] then
            Out[i].n_flow[j] = beta[i, j] * o[j];
          end if;
        end for;
      end for;
      x[:, 1] = Out[2].n_flow ./ sum(Out[2].n_flow);
      y[:, 1] = Out[1].n_flow ./ sum(Out[1].n_flow);
      if known_P then
        P = Pressure;
      else
        for i in 1:nin loop
          dummy_sum[i] = In[i].P;
        end for;
        P = sum(dummy_sum) / nin;
      end if;
      if system.use_T then
        if known_T then
          T = Temperature;
          phiL = PR_EOS(system.set_n, x, P, T, system.pc, system.Tc, system.w, system.kij, "L");
          phiV = PR_EOS(system.set_n, y, P, T, system.pc, system.Tc, system.w, system.kij, "V");
          K = phiL ./ phiV;
          y = x .* K;
        else
          (T, K) = PTFlash(system.set_n, x, y, P, T_guess, dT, system.pc, system.Tc, system.w, system.kij);
          y = x .* K;
        end if;
        for i in 1:nout loop
          Out[i].T = T;
        end for;
      end if;
      annotation(Icon(coordinateSystem(extent = {{-55, -125}, {100, 125}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Line(points = {{-50, 60}, {50, -60}}, smooth = Smooth.Bezier), Line(points = {{-50, -60}, {50, 60}}, smooth = Smooth.Bezier), Line(origin = {38.993, 52.3419}, points = {{-38.993, 37.5878}, {-38.993, 52.3419}, {38.993, 52.3419}, {38.993, -6.67447}, {10.1874, -6.67447}, {10.8899, -6.67447}}), Line(origin = {38.993, -52.693}, points = {{-38.993, -37.9393}, {-38.993, -50.5857}, {38.993, -50.5857}, {38.993, 7.72815}, {10.8899, 7.72815}, {10.8899, 7.72815}}), Line(origin = {78.88, -103.98}, points = {{-21.7798, 21.0773}, {19.6721, -20.3747}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}, smooth = Smooth.Bezier), Line(origin = {78.2, 105.58}, rotation = 90, points = {{-23.185, 22.4824}, {19.6721, -20.3747}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}, smooth = Smooth.Bezier), Rectangle(origin = {76, 103}, fillColor = {198, 200, 203}, fillPattern = FillPattern.VerticalCylinder, extent = {{-15, -15}, {15, 15}}, radius = 15), Rectangle(origin = {77.02, -102.32}, fillColor = {198, 200, 203}, fillPattern = FillPattern.VerticalCylinder, extent = {{-15, -15}, {15, 15}}, radius = 15)}));
    end Sep2;

    model RStoic
      extends Chemical.Icons.Column;
      import Chemical.Utilities.*;
      parameter Integer nin = 1 "Number of inputs" annotation(Dialog(group = "Number of Ports"));
      parameter Integer nout = 2 "Number of Output" annotation(Dialog(group = "Number of Ports", enable = false));
      extends Chemical.Reactive_Sys.Reactor_generic(final Bernoulli = false);
      parameter Integer[n_rxn] BC = {1, 2, 1} "Enter the component index of Base Component foe all the reactions" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter Real[n_rxn] f = 0 "Enter fractional conversion of BC for all the reactions, for eg. {0.9,0.8,0.2}" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter Boolean Selectivity = false "Enter true if calculations are to be based on Selectivity" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter Real Sp = 1 "Enter the selectivity of desired product w.r.t. Base Component" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter Integer nY = 2 "Enter the component index of desired product" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter Boolean known_T "Give true if operating Temperature is known" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter Real Temperature "Enter operating Temperature" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter Real T_guess = 273 "Enter guess valur for equilibrium Temperature if operating temperature is not known" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      Real Y_max if Selectivity "Maximum possible desired product";
      Real Y if Selectivity;
      Real[n_rxn] flag1 if Selectivity;
      Real P;
      Real T if system.use_T;
      constant Real dT = 20;
      constant Real dT = 20;
      Real[system.set_n, 1] y(start = ones(system.set_n, 1)) if not SinglePhase;
      Real[system.set_n, 1] x(start = ones(system.set_n, 1)) if not SinglePhase;
      Real[system.set_n, 1] K if not SinglePhase;
      parameter Boolean SinglePhase "Give true if the system has only single phase" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter String ph = "L" "Enter phase of Reaction, L for Liquid and V for Gas" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter Real[system.set_n, 1] output_guess = ones(system.set_n, 1) "Enter guess values for the outlet mixture" annotation(Dialog(tab = "Reactor Spec", group = "Product"));
      Chemical.Interfaces.Mports_In In[nin](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {-99.5, -2}, extent = {{-12.5, -50}, {12.5, 50}}, rotation = 0), iconTransformation(origin = {-50, 0}, extent = {{-5, -20}, {5, 20}}, rotation = 0)));
      Chemical.Interfaces.Columnport Out[nout](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {0, 0}, extent = {{-23.2812, -93.125}, {23.2812, 93.125}}, rotation = 0), iconTransformation(origin = {0.5, 3.55271e-015}, extent = {{-13.5, -90}, {13.5, 90}}, rotation = 0)));
    equation
      if Selectivity then
        Y = Sp;
        if not coeff[1, BC] == 0 then
          flag1[1] = coeff[1, nY] / (-coeff[1, BC]);
        else
          flag1[1] = 1;
        end if;
        for i in 2:n_rxn loop
          if not coeff[i, BC] == 0 then
            flag1[i] = flag[i - 1] * coeff[i, nY] / (-coeff[i, BC]);
          else
            flag1[i] = 1;
          end if;
        end for;
        Y_max = o[BC] * flag1[n_rxn];
        p[nY] = Y * Y_max;
        for i in 1:n_rxn - 1 loop
          if not coeff[n_rxn, nY] == 0 then
            e[i, 1] = (p[nY] - o[nY]) / coeff[n_rxn, nY];
          end if;
        end for;
      else
        for i in 1:n_rxn loop
          e[i, 1] = -f[i] * o[BC[i]] / coeff[i, BC[i]];
        end for;
      end if;
      // Flash Calculations
      P = In[1].P - Pdrop;
      Out[1].P = P;
      Out[2].P = P;
      if system.use_T then
        if not SinglePhase then
          x[:, 1] = Out[2].n_flow ./ sum(Out[2].n_flow);
          y[:, 1] = Out[1].n_flow ./ sum(Out[1].n_flow);
          if known_T then
            T = Temperature;
            phiL1 = PR_EOS(system.set_n, x, P, T, system.pc, system.Tc, system.w, system.kij, "L");
            phiV1 = PR_EOS(system.set_n, y, P, T, system.pc, system.Tc, system.w, system.kij, "V");
            K = phiL1 ./ phiV1;
            y = x .* K;
          else
            (T, K) = PTFlash(system.set_n, x, y, P, T_guess, dT, system.pc, system.Tc, system.w, system.kij, ph);
            y = x .* K;
          end if;
        else
          if known_T then
            T = Temperature;
          else
            T = Out[1].T;
          end if;
          if ph == "L" then
            for i in 1:system.set_n loop
              Out[1].n_flow[i] = 0;
            end for;
          elseif ph == "V" then
            for i in 1:system.set_n loop
              Out[2].n_flow[i] = 0;
            end for;
          end if;
        end if;
        Out[1].T = Out[2].T;
      end if;
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-55, -90}, {50, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Line(origin = {1, 19}, points = {{-50.9541, -0.318021}, {50, -0.318021}}), Line(origin = {1, -19}, points = {{-50, 0}, {50, 0}}), Line(points = {{-50, -19}, {-19, 19}}, thickness = 1, smooth = Smooth.Bezier), Line(points = {{-50, 19}, {-19, -19}}, thickness = 1, smooth = Smooth.Bezier), Line(origin = {31.4488, -0.0353357}, points = {{-50, -19}, {-19, 19}}, thickness = 1, smooth = Smooth.Bezier), Line(origin = {30.8127, -0.0353357}, points = {{-50, 19}, {-19, -19}}, thickness = 1, smooth = Smooth.Bezier), Line(origin = {61.9435, -0.388693}, points = {{-50, -19}, {-19, 19}}, thickness = 1, smooth = Smooth.Bezier), Line(origin = {62.2615, 0.24735}, points = {{-50, 19}, {-19, -19}}, thickness = 1, smooth = Smooth.Bezier)}));
    end RStoic;

    model REquil
      extends Chemical.Icons.Column;
      import Chemical.Utilities.*;
      parameter Integer nin = 1 "Number of inputs" annotation(Dialog(group = "Number of Ports"));
      parameter Integer nout = 2 "Number of Output" annotation(Dialog(group = "Number of Ports", enable = false));
      extends Chemical.Reactive_Sys.Reactor_generic(final Bernoulli = false);
      Chemical.Interfaces.Mports_In In[nin](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {-99.5, -2}, extent = {{-12.5, -50}, {12.5, 50}}, rotation = 0), iconTransformation(origin = {-50, 0}, extent = {{-5, -20}, {5, 20}}, rotation = 0)));
      parameter Boolean known_T "Give true if operating Temperature is known" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter Real Temperature "Enter operating Temperature" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter Real T_guess = 273 "Enter guess valur for equilibrium Temperature if operating temperature is not known" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter Real[n_rxn] Kequ "Enter Equillibrium constant for all the reactions" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter Boolean SinglePhase "Give true if the system has only single phase" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter String ph = "L" "Enter phase of Reaction, L for Liquid and V for Gas" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter Real[system.set_n, 1] output_guess = ones(system.set_n, 1) "Enter guess values for the outlet mixture" annotation(Dialog(tab = "Reactor Spec", group = "Product"));
      Real[n_rxn, system.set_n] dummy_product "for calculation only. Please ignore";
      Real P;
      Real T if system.use_T;
      constant Real dT = 20;
      Real[system.set_n, 1] y(start = ones(system.set_n, 1)) if not SinglePhase;
      Real[system.set_n, 1] x(start = ones(system.set_n, 1)) if not SinglePhase;
      Real[system.set_n, 1] z(start = output_guess) if system.use_T;
      Real[system.set_n, 1] phiL1 if not SinglePhase and known_T;
      Real[system.set_n, 1] phiV1 if not SinglePhase and known_T;
      Real[system.set_n, 1] phi if system.use_T;
      Real[system.set_n, 1] K(start = ones(system.set_n, 1)) if not SinglePhase;
      Chemical.Interfaces.Columnport Out[nout](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {0, 0}, extent = {{-23.2812, -93.125}, {23.2812, 93.125}}, rotation = 0), iconTransformation(origin = {0.5, 3.55271e-015}, extent = {{-13.5, -90}, {13.5, 90}}, rotation = 0)));
    equation
      P = In[1].P - Pdrop;
      Out[1].P = P;
      Out[2].P = P;
      if system.use_T then
        if not SinglePhase then
          x[:, 1] = Out[2].n_flow ./ sum(Out[2].n_flow);
          y[:, 1] = Out[1].n_flow ./ sum(Out[1].n_flow);
          if known_T then
            T = Temperature;
            phiL1 = PR_EOS(system.set_n, x, P, T, system.pc, system.Tc, system.w, system.kij, "L");
            phiV1 = PR_EOS(system.set_n, y, P, T, system.pc, system.Tc, system.w, system.kij, "V");
            K = phiL1 ./ phiV1;
            y = x .* K;
          else
            (T, K) = PTFlash(system.set_n, x, y, P, T_guess, dT, system.pc, system.Tc, system.w, system.kij, ph);
            y = x .* K;
          end if;
        else
          if known_T then
            T = Temperature;
          else
            T = Out[1].T;
          end if;
          if ph == "L" then
            for i in 1:system.set_n loop
              Out[1].n_flow[i] = 0;
            end for;
          elseif ph == "V" then
            for i in 1:system.set_n loop
              Out[2].n_flow[i] = 0;
            end for;
          end if;
        end if;
        Out[1].T = Out[2].T;
        z[:, 1] = p ./ sum(p);
        phi = PR_EOS(system.set_n, z, P, T, system.pc, system.Tc, system.w, system.kij, ph);
      end if;
      for i in 1:n_rxn loop
        dummy_product[i, 1] = (phi[1, 1] * z[1, 1] * P / 101325) ^ coeff[i, 1];
        for j in 2:system.set_n loop
          dummy_product[i, j] = dummy_product[i, j - 1] * (phi[j, 1] * z[j, 1] * P / 101325) ^ coeff[i, j];
        end for;
      end for;
      Kequ = dummy_product[:, system.set_n];
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-55, -90}, {50, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Line(points = {{0, -60}, {0, 60}}, thickness = 0.75), Line(origin = {-20, 0}, points = {{0, -60}, {0, 60}}, thickness = 0.75), Line(origin = {10, 0}, points = {{0, -60}, {0, 60}}, thickness = 0.75), Line(origin = {-10, 0}, points = {{0, -60}, {0, 60}}, thickness = 0.75), Line(origin = {-30, 0}, points = {{0, -60}, {0, 60}}, thickness = 0.75), Line(origin = {20, 0}, points = {{0, -60}, {0, 60}}, thickness = 0.75), Line(origin = {30, 0}, points = {{0, -60}, {0, 60}}, thickness = 0.75), Line(origin = {40, 0}, points = {{0, -60}, {0, 60}}, thickness = 0.75), Line(origin = {-40, 0}, points = {{0, -60}, {0, 60}}, thickness = 0.75)}));
    end REquil;

    model Cooler
      extends Icons.Heater;
      import Chemical.Utilities.*;
      parameter Integer nin = 1 "Number of input(s)" annotation(Dialog(group = "Number of Ports", enable = false));
      parameter Integer nout = 1 "Number of output(s)" annotation(Dialog(group = "Number of Ports", enable = false));
      extends Chemical.NonReactive_Sys.GeneralBal(final known_W = true, final Work = 0, final Bernoulli = false, final known_FLoss = false, final F_loss = 0, final Temp_change = true);
      parameter Boolean known_P = true "Give true if Outlet Presuure is known" annotation(Dialog(tab = "Mode"));
      parameter Real P "Enter the Outlet Pressure, if known" annotation(Dialog(tab = "Mode"));
      parameter Boolean known_Tout = false "Give true if Outlet Temperature is known" annotation(Dialog(tab = "Mode"));
      parameter Real Tout "Enter the Outlet Temperature, if known" annotation(Dialog(tab = "Mode"));
      parameter Boolean known_Vfraction = false "Give true if Vapour Phase molar fraction is known" annotation(Dialog(tab = "Mode"));
      parameter Real Vfraction "Enter Vapour Phase molar fraction" annotation(Dialog(tab = "Mode"));
      parameter Boolean Heat_Stream = false "Give true if Heat stream is used" annotation(Dialog(tab = "Mode"));
      Chemical.Interfaces.port_In In[nin](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {-80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Chemical.Interfaces.port_Out Out[nout](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {79, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Chemical.Interfaces.HeatPort_In QIn if Heat_Stream annotation(Placement(visible = true, transformation(origin = {0, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      if known_P then
        Out[1].P = P;
      end if;
      if not known_Q then
        if known_Tout then
          Out[1].T = Tout;
        elseif known_Vfraction then
          sum(Out[1].n_flow) = Vfraction * sum(o);
        end if;
      end if;
      if Heat_Stream then
        QIn.Q_flow = Q;
      end if;
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Text(lineColor = {77, 99, 179}, extent = {{-60, 60}, {60, -60}}, textString = "C", textStyle = {TextStyle.Bold})}));
    end Cooler;

    model Heater
      extends Icons.Heater;
      import Chemical.Utilities.*;
      parameter Integer nin = 1 "Number of input(s)" annotation(Dialog(group = "Number of Ports", enable = false));
      parameter Integer nout = 1 "Number of output(s)" annotation(Dialog(group = "Number of Ports", enable = false));
      extends Chemical.NonReactive_Sys.GeneralBal(final Bernoulli = false, final known_W = true, final Work = 0, final known_FLoss = false, final F_loss = 0, final Temp_change = true);
      parameter Boolean known_P = true "Give true if Outlet Presuure is known" annotation(Dialog(tab = "Mode"));
      parameter Real P "Enter the Outlet Pressure, if known" annotation(Dialog(tab = "Mode"));
      parameter Boolean known_Tout = false "Give true if Outlet Temperature is known" annotation(Dialog(tab = "Mode"));
      parameter Real Tout "Enter the Outlet Temperature, if known" annotation(Dialog(tab = "Mode"));
      parameter Boolean known_Vfraction = false "Give true if Vapour Phase molar fraction is known" annotation(Dialog(tab = "Mode"));
      parameter Real Vfraction "Enter Vapour Phase molar fraction" annotation(Dialog(tab = "Mode"));
      parameter Boolean Heat_Stream = false "Give true if Heat stream is used" annotation(Dialog(tab = "Mode"));
      Chemical.Interfaces.port_In In[nin](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {-80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Chemical.Interfaces.port_Out Out[nout](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {79, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Chemical.Interfaces.HeatPort_Out QOut if Heat_Stream annotation(Placement(visible = true, transformation(origin = {0, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      if known_P then
        Out[1].P = P;
      end if;
      if not known_Q then
        if known_Tout then
          Out[1].T = Tout;
        elseif known_Vfraction then
          sum(Out[1].n_flow) = Vfraction * sum(o);
        end if;
      end if;
      if Heat_Stream then
        QOut.Q_flow = Q;
      end if;
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Text(origin = {0, -0.468384}, lineColor = {220, 0, 0}, extent = {{-60, 60}, {60, -60}}, textString = "H", textStyle = {TextStyle.Bold})}));
    end Heater;

    model ShortcutHE
      parameter Boolean known_Qc = false "Give true if the heat transferred to/from the system is known" annotation(Dialog(tab = "Cooler", group = "Heat"));
      parameter Real Heatflow_C "Give Heat transf. to the system as +ve else if the system is adiabatic = 0" annotation(Dialog(tab = "Cooler", group = "Heat"));
      parameter Boolean known_Tout_C = false "Give true if Outlet Temperature is known" annotation(Dialog(tab = "Cooler", group = "Mode"));
      parameter Real Tout_C "Enter the Outlet Temperature, if known" annotation(Dialog(tab = "Cooler", group = "Mode"));
      parameter Boolean known_VfractionC = false "Give true if Vapour Phase molar fraction is known" annotation(Dialog(tab = "Cooler", group = "Mode"));
      parameter Real Vfraction_C "Enter Vapour Phase molar fraction" annotation(Dialog(tab = "Cooler", group = "Mode"));
      parameter Boolean known_Qh = false "Give true if the heat transferred to/from the system is known" annotation(Dialog(tab = "Heater", group = "Heat"));
      parameter Real Heatflow_H "Give Heat transf. to the system as +ve else if the system is adiabatic = 0" annotation(Dialog(tab = "Heater", group = "Heat"));
      parameter Boolean known_Tout_H = false "Give true if Outlet Temperature is known" annotation(Dialog(tab = "Heater", group = "Mode"));
      parameter Real Tout_H "Enter the Outlet Temperature, if known" annotation(Dialog(tab = "Heater", group = "Mode"));
      parameter Boolean known_VfractionH = false "Give true if Vapour Phase molar fraction is known" annotation(Dialog(tab = "Heater", group = "Mode"));
      parameter Real Vfraction_H "Enter Vapour Phase molar fraction" annotation(Dialog(tab = "Heater", group = "Mode"));
      Chemical.UnitOp.Cooler C1(final Heat_Stream = true, known_Q = known_Qc, Heatflow = Heatflow_C, known_Tout = known_Tout_C, Tout = Tout_C, known_Vfraction = known_VfractionC, Vfraction = Vfraction_C) annotation(Placement(visible = true, transformation(origin = {0, 60}, extent = {{-43.75, -43.75}, {43.75, 43.75}}, rotation = 0)));
      Chemical.UnitOp.Heater H1(final Heat_Stream = true, known_Q = known_Qh, Heatflow = Heatflow_H, known_Tout = known_Tout_H, Tout = Tout_H, known_Vfraction = known_VfractionH, Vfraction = Vfraction_H) annotation(Placement(visible = true, transformation(origin = {1.5, -50.5}, extent = {{47.5, -47.5}, {-47.5, 47.5}}, rotation = 0)));
      Chemical.Interfaces.port_In port_inH annotation(Placement(visible = true, transformation(origin = {100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Chemical.Interfaces.port_Out port_outH annotation(Placement(visible = true, transformation(origin = {-100, -49}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -47}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Chemical.Interfaces.port_In portInC annotation(Placement(visible = true, transformation(origin = {-100, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Chemical.Interfaces.port_Out portOutC annotation(Placement(visible = true, transformation(origin = {101, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(portInC, C1.In) annotation(Line(points = {{-100, 60}, {-34.6604, 60}, {-34.6604, 60}, {-35, 60}}, color = {0, 0, 255}, arrow = {Arrow.None, Arrow.Filled}));
      connect(port_outH, H1.Out) annotation(Line(points = {{-100, -49}, {-37.4707, -49}, {-37.4707, -48.6}, {-36.025, -48.6}}, color = {0, 0, 255}));
      connect(H1.In, port_inH) annotation(Line(points = {{39.5, -50.5}, {100.234, -50.5}, {100.234, -50}, {100, -50}}, color = {0, 0, 255}));
      connect(C1.Out, portOutC) annotation(Line(points = {{34.5625, 61.75}, {101.639, 61.75}, {101.639, 62}, {101, 62}}, color = {0, 0, 255}));
      connect(heater1.QOut, cooler1.QIn) annotation(Line(points = {{1.5, -12.5}, {0.468384, -12.5}, {0.468384, 25}, {0, 25}}, pattern = LinePattern.Dash, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 5));
      annotation(Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(extent = {{-100, -26}, {100, -30}}, lineColor = {0, 0, 0}, fillColor = {95, 95, 95}, fillPattern = FillPattern.Forward), Rectangle(extent = {{-100, 30}, {100, 26}}, lineColor = {0, 0, 0}, fillColor = {95, 95, 95}, fillPattern = FillPattern.Forward), Rectangle(extent = {{-100, 60}, {100, 30}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.HorizontalCylinder, fillColor = {0, 63, 125}), Rectangle(extent = {{-100, -30}, {100, -60}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.HorizontalCylinder, fillColor = {0, 63, 125}), Rectangle(extent = {{-100, 26}, {100, -26}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.HorizontalCylinder, fillColor = {0, 128, 255}), Text(extent = {{-150, 110}, {150, 70}}, lineColor = {0, 0, 255}, textString = ""), Line(points = {{30, -85}, {-60, -85}}, color = {0, 128, 255}, smooth = Smooth.None), Polygon(points = {{20, -70}, {60, -85}, {20, -100}, {20, -70}}, lineColor = {0, 128, 255}, smooth = Smooth.None, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Line(points = {{30, 77}, {-60, 77}}, color = {0, 128, 255}, smooth = Smooth.None), Polygon(points = {{-50, 92}, {-90, 77}, {-50, 62}, {-50, 92}}, lineColor = {0, 128, 255}, smooth = Smooth.None, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid)}));
    end ShortcutHE;

    model DSTFUG
      extends Icons.Column;
      import Chemical.Utilities.*;
      extends Chemical.NonReactive_Sys.GeneralBal(final known_W = true, final Work = 0, final Bernoulli = false, final known_FLoss = false, final F_loss = 0, final known_Q = false, final Heatflow = 0, final Temp_change = true);
      parameter Integer nin = 1 "Number of inputs" annotation(Dialog(group = "Number of Ports", enable = false));
      parameter Integer nout = 2 "Number of Output" annotation(Dialog(group = "Number of Ports", enable = false));
      Chemical.Interfaces.port_In In[nin](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {-60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-49, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Chemical.Interfaces.Columnport Out[nout](n = system.set_n, EnergyBal = system.use_EnergyBal, OpenSys = system.OpenSys, use_T = system.use_T) annotation(Placement(visible = true, transformation(origin = {0, 0}, extent = {{-23.2812, -93.125}, {23.2812, 93.125}}, rotation = 0), iconTransformation(origin = {0.5, 3.55271e-015}, extent = {{-13.5, -90}, {13.5, 90}}, rotation = 0)));
      parameter Boolean known_N "Give true if number of theoretical stages are known" annotation(Dialog(tab = "Specifications", group = "Column"));
      parameter Real N_theo "Enter number of theoretical stages for the column" annotation(Dialog(tab = "Specifications", group = "Column"));
      parameter Boolean known_R "Give true if molar refux ratio is known" annotation(Dialog(tab = "Specifications", group = "Column"));
      parameter Real RefluxRatio "Enter molar reflux ratio" annotation(Dialog(tab = "Specifications", group = "Column"));
      parameter Real CondenserP "Enter Operating pressure of the condenser" annotation(Dialog(tab = "Specifications", group = "Pressure"));
      parameter Real ReboilerP "Enter Operating pressure of the reboiler" annotation(Dialog(tab = "Specifications", group = "Pressure"));
      parameter Integer Index_LK "Enter index of the Light key component" annotation(Dialog(tab = "Specifications", group = "Component Recoveries"));
      parameter Real xLK "Enter Light key molar fraction in Bottoms" annotation(Dialog(tab = "Specifications", group = "Component Recoveries"));
      parameter Integer Index_HK "Enter index of the Heavy key component" annotation(Dialog(tab = "Specifications", group = "Component Recoveries"));
      parameter Real xHK "Enter Heavy key molar fraction in Distillate" annotation(Dialog(tab = "Specifications", group = "Component Recoveries"));
      parameter Integer CType(min = 1, max = 2) "Give 1 if it is a total condenser, 2 for a partial condenser with all vapour distillate" annotation(Dialog(tab = "Specifications", group = "Condenser"));
      Real N "Theoretical number of stages";
      Real Nmin "Minimum number of theoretical stages";
      Real R "Actual Reflux ratio";
      Real Rmin "Minimum Reflux ratio";
      Real Td "Distillate Temperature";
      Real Tb "Bottoms Temperature";
      Real[system.set_n, 1] Kd;
      Real[system.set_n, 1] Kb;
      Real[system.set_n, 1] Kbubl;
      Real theta;
      Real q "Heat Required to vaporize 1 mole of feed/Avg. Latent Heat of Feed";
      Real Tbubl "Bubble point temperature of the feed";
      Real[system.set_n, 4] H1;
      Real Q1;
      Real Qc "Condenser Duty";
      Real Qr "Reboiler Duty";
      Real[system.set_n, 1] Tr "Reduced Temperature";
      Real[system.set_n, 1] Tr_d "Reduced Temperature for distillate";
      Real[system.set_n, 1] Tr_b "Reduced Temperature for bottoms";
      Real[system.set_n, 1] lambda;
      Real[system.set_n, 1] lambda_d;
      Real[system.set_n, 1] lambda_b;
      Real[systems.set_n, 1] alphad;
      Real[systems.set_n, 1] alphab;
      constant Real dT = 20;
      Real[system.set_n, 1] alpha;
      Real[system.set_n] dummy;
      Real[system.set_n] dummy2;
      Real dummy3;
      Real dummy4;
      Real[system.set_n, 1] x;
      Real[system.set_n, 1] y;
      Real[system.set_n, 1] z;
    equation
      x[:, 1] = Out[2].n_flow ./ sum(Out[2].n_flow);
      y[:, 1] = Out[1].n_flow ./ sum(Out[1].n_flow);
      z[:, 1] = In[1].n_flow ./ sum(In[1].n_flow);
      (Tbubl, Kbubl) = VLE_PR(system.set_n, z, In[1].P, In[1].T, dT, system.pc, system.Tc, system.w, system.kij);
      Tr = In[1].T .* ones(system.set_n, 1) ./ system.Tc;
      for i in 1:system.set_n loop
        H1[i, 1] = system.Cp_coeff_gas[i, 1] * (Tbubl - port.T);
        lambda[i, 1] = system.delHvap[i, 1] * (1 - Tr[i, 1]) ^ (delHvap[i, 2] + delHvap[i, 3] * Tr[i, 1] + delHvap[i, 4] * Tr[i, 1] ^ 2);
        for j in 2:4 loop
          H1[i, j] = H1[i, j - 1] + system.Cp_coeff_gas[i, j] / j * (Tbubl ^ j - port.T ^ j);
        end for;
      end for;
      Q1 = sum((transpose(H1[:, 4]) + lambda) .* o / sum(o));
      q = Q1 / sum(lambda .* o / sum(o));
      if known_N then
        N = N_theo;
      elseif known_R then
        R = RefluxRatio;
      end if;
      Out[1].n_flow[Index_HK] / sum(Out[1].n_flow) = xHK;
      Out[2].n_flow[Index_LK] / sum(Out[2].n_flow) = xLK;
      (Td, Kd) = VLE_PR(system.set_n, y, CondenserP, In[1].T, dT, system.pc, system.Tc, system.w, system.k);
      (Tb, Kb) = VLE_PR(system.set_n, x, ReboilerP, In[1].T, dT, system.pc, system.Tc, system.w, system.k);
      Out[1].T = Td;
      Out[2].T = Tb;
      alphad = Kd ./ Kd[Index_HK, 1];
      alphab = Kb ./ Kb[Index_HK, 1];
      alpha = sqrt(alphad .* alphab);
      // FENSKE EQUATION
      for i in 1:system.set_n loop
        if not i == Index_HK then
          Nmin = log(y[i] * x[Index_HK] / (y[Index_HK] * x[i])) / log(alpha[i, 1]);
        end if;
      end for;
      // UNDERWOOD CORRELATION
      for i in 1:system.set_n loop
        dummy[i] = alpha[i, 1] * z[i] / (alpha[i, 1] - theta);
        dummy2[i] = alpha[i, 1] * y[i] / (alpha[i, 1] - theta);
      end for;
      1 - q = sum(dummy);
      Rmin + 1 = sum(dummy2);
      // GILLILAND CORRELATION
      dummy3 = R - Rmin / (R + 1);
      dummy4 = N - Nmin / (N + 1);
      dummy4 = 0.75 * (1 - dummy3 ^ 0.5668);
      // Heat Duties
      Q = Qc + Qr;
      Tr_d = Td .* ones(system.set_n, 1) ./ system.Tc;
      Tr_b = Tb .* ones(system.set_n, 1) ./ system.Tc;
      for i in 1:system.set_n loop
        lambda_d[i, 1] = system.delHvap[i, 1] * (1 - Tr_d[i, 1]) ^ (delHvap[i, 2] + delHvap[i, 3] * Tr_d[i, 1] + delHvap[i, 4] * Tr_d[i, 1] ^ 2);
        lambda_b[i, 1] = system.delHvap[i, 1] * (1 - Tr_b[i, 1]) ^ (delHvap[i, 2] + delHvap[i, 3] * Tr_b[i, 1] + delHvap[i, 4] * Tr_b[i, 1] ^ 2);
      end for;
      if CType == 1 then
        Qc = sum(lambda_d .* Out[1].n_flow);
      elseif CType == 2 then
        Qc = sum(lambda_d .* R .* Out[1].n_flow);
      end if;
      Qr = sum(lambda_b .* Out[2].n_flow);
      annotation(Icon(coordinateSystem(extent = {{-55, -125}, {100, 125}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Line(origin = {38.99, 52.34}, points = {{-38.993, 37.5878}, {-38.993, 52.3419}, {38.993, 52.3419}, {38.993, -6.67447}, {10.1874, -6.67447}, {10.8899, -6.67447}}), Line(origin = {38.993, -52.693}, points = {{-38.993, -37.9393}, {-38.993, -50.5857}, {38.993, -50.5857}, {38.993, 7.72815}, {10.8899, 7.72815}, {10.8899, 7.72815}}), Line(origin = {78.88, -103.98}, points = {{-21.7798, 21.0773}, {19.6721, -20.3747}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}, smooth = Smooth.Bezier), Line(origin = {78.2, 105.58}, rotation = 90, points = {{-23.185, 22.4824}, {19.6721, -20.3747}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}, smooth = Smooth.Bezier), Rectangle(origin = {76, 103}, fillColor = {198, 200, 203}, fillPattern = FillPattern.VerticalCylinder, extent = {{-15, -15}, {15, 15}}, radius = 15), Rectangle(origin = {77.02, -102.32}, fillColor = {198, 200, 203}, fillPattern = FillPattern.VerticalCylinder, extent = {{-15, -15}, {15, 15}}, radius = 15), Line(origin = {-0.04, 37.81}, points = {{-49.47, 4}, {50, 4}}), Line(origin = {-0.181343, 20.5273}, points = {{-49.47, 4}, {50, 4}}), Line(origin = {0.20735, 3.18456}, points = {{-49.47, 4}, {50, 4}}), Line(origin = {0.0059364, -12.6282}, points = {{-49.47, 4}, {50, 4}}), Line(origin = {0.5, -27.33}, points = {{-50, 0}, {48.9399, 0}}), Line(origin = {0.358657, -47.1427}, points = {{-50, 0}, {48.9399, 0}})}));
    end DSTFUG;
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end UnitOp;

  package Icons
    extends Modelica.Icons.IconsPackage;

    model Column
      annotation(Icon(coordinateSystem(extent = {{-50, -90}, {50, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(fillColor = {230, 225, 231}, fillPattern = FillPattern.VerticalCylinder, extent = {{-50, 60}, {50, -60}}), Ellipse(origin = {0, 70}, fillColor = {198, 200, 203}, fillPattern = FillPattern.Solid, extent = {{-50, -40}, {50, 20}}, endAngle = 180), Ellipse(origin = {0, -50}, fillColor = {198, 200, 203}, fillPattern = FillPattern.Solid, extent = {{-50, -40}, {50, 20}}, endAngle = -180)}));
    end Column;

    model Mixer
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Polygon(fillColor = {216, 222, 221}, fillPattern = FillPattern.HorizontalCylinder, points = {{-100, 100}, {60, 20}, {100, 20}, {100, -20}, {60, -20}, {-100, -100}, {-100, 100}})}));
    end Mixer;

    model Splitter
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Polygon(rotation = 180, fillColor = {216, 222, 221}, fillPattern = FillPattern.VerticalCylinder, points = {{-100, 100}, {60, 20}, {100, 20}, {100, -20}, {60, -20}, {-100, -100}, {-100, 100}})}));
    end Splitter;

    model Heater
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-80, -80}, {80, 80}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Ellipse(fillColor = {198, 200, 203}, fillPattern = FillPattern.Solid, extent = {{-80, 80}, {80, -80}}, endAngle = 360), Line(points = {{-60, 0}, {-30, 30}, {30, -30}, {60, 0}}), Line(origin = {2, 0}, points = {{-80, -80}, {80, 80}}, thickness = 0.75, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 10, smooth = Smooth.Bezier), Line(points = {{-60, 0}, {-80, 0}}), Line(points = {{60, 0}, {80, 0}})}));
    end Heater;
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end Icons;

  package Reactive_Sys
    extends Modelica.Icons.VariantsPackage;

    model Reactor_generic
      import Chemical.Utilities.*;
      extends Chemical.GenericModel;
      parameter Boolean Bernoulli = false "Give true if Heat changes are negligible" annotation(Dialog(tab = "Energy Transfer", group = "Mechanical"));
      parameter Integer n_rxn(start = 1) "number of reactions" annotation(Dialog(tab = "Reactor Spec", group = "Stoichiometry"));
      parameter Real[n_rxn, system.set_n] coeff "Stoichiometric coefficients of the components in the different reactions. e.g. {{-1,1,0},{-1,0,2}}" annotation(Dialog(tab = "Reactor Spec", group = "Stoichiometry"));
      parameter Real[n_rxn] H_rxn "Enter the enthalpies of reaction" annotation(Dialog(tab = "Reactor Spec", group = "Enthalpy"));
      parameter Boolean known_Q "Give true if the heat transferred to/from the system is known" annotation(Dialog(tab = "Energy Transfer", group = "Heat"));
      parameter Real Heatflow "Give Heat transf. to the system as +ve else if the system is adiabatic = 0" annotation(Dialog(tab = "Energy Transfer", group = "Heat"));
      parameter Boolean Isothermic "Give true if the  system is Isothermic" annotation(Dialog(tab = "Energy Transfer", group = "Heat"));
      parameter Real Pdrop = 0 "Enter Pressure drop across the reactor" annotation(Dialog(tab = "Reactor Spec", group = "Reaction Parameters"));
      parameter Real[n_rxn, 1] extent_guess = ones(n_rxn, 1) "Enter guess values for the outlet mixture" annotation(Dialog(tab = "Reactor Spec", group = "Product"));
      Real Q if system.use_EnergyBal and not Bernoulli "Heat transferred to/from the system from/to the surroundings";
      Real[n_rxn] delHr;
      //(start = {{0.34326}, {1.54865}})
      Real[n_rxn, 1] e(start = extent_guess) "Extent of Reaction";
      Real[n_rxn, system.set_n] dummy_sum "for calculation only. Please ignore";
      Integer LR(start = 1) "Limiting Reagent";
    equation
      LR = Limiting_Reagent(system.set_n, o, coeff[1, :]);
      dummy_sum[1, :] = coeff[1, :] * e[1, 1];
      for j in 2:n_rxn loop
        dummy_sum[j, :] = dummy_sum[j - 1, :] + coeff[j, :] * e[j, 1];
      end for;
      p = o + dummy_sum[n_rxn, :];
      // Energy Balance
      if system.use_EnergyBal then
        if Isothermic then
          T = In[1].T;
        end if;
        if known_Q then
          Q = Heatflow;
        end if;
        for i in 1:n_rxn loop
          delHr[i] = e[i, 1] * H_rxn[i];
        end for;
        delH + delEk + delEp = Q + sum(delHr);
      end if;
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
    end Reactor_generic;

    model R_Combustion
      extends Reactor_Stoic;
      parameter Boolean known_xsair(start = false) "Give true if percentage excess air is known";
      parameter Real xsair(start = 0) "Enter percentage excess air";
      parameter Integer O2_Ind "Give component index of O2";
      Real xsair1 "Only for calculation";
      Real air_th "Theoretical air";
    equation
      if known_xsair then
        xsair1 = xsair;
      end if;
      air_th = In[1].n_flow[1] * (-coeff[1, O2_Ind]);
      In[2].n_flow[O2_Ind] = (1 + xsair1 / 100) * air_th;
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(origin = {0.936768, -3.55271e-015}, fillColor = {241, 15, 30}, fillPattern = FillPattern.VerticalCylinder, extent = {{-100, 100}, {100, -100}})}));
    end R_Combustion;
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end Reactive_Sys;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
end Chemical;