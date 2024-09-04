// Copyright 2023 ETH Zurich and University of Bologna.
// Solderpad Hardware License, Version 0.51, see LICENSE for details.
// SPDX-License-Identifier: SHL-0.51
//
// Description: Automatically generated bootrom
//
// Generated by hardware/scripts/generate_bootrom.py

module bootrom #(
  /* Automatically generated. DO NOT CHANGE! */
  parameter int unsigned DataWidth = 64,
  parameter int unsigned AddrWidth = 48
) (
  input  logic                 clk_i,
  input  logic                 req_i,
  input  logic [AddrWidth-1:0] addr_i,
  output logic [DataWidth-1:0] rdata_o
);
  localparam int RomSize = 17;
  localparam int AddrBits = RomSize > 1 ? $clog2(RomSize) : 1;

  const logic [RomSize-1:0][DataWidth-1:0] mem = {
    64'h00000000ffffffff,
    64'h0000104800001044,
    64'h0000000000001044,
    64'h0000000078400000,
    64'h0000000078000000,
    64'h0000000000020000,
    64'h5100000000000010,
    64'h0000000200001000,
    64'h0000006f00038067,
    64'h0003a38300038393,
    64'h0583839301c383b3,
    64'h0105ae0300c5a383,
    64'h1050007330461073,
    64'h0086661330402673,
    64'h06c5a58300000597,
    64'hf140257330531073,
    64'h0783230300000317
  };

  logic [AddrBits-1:0] addr_q;

  always_ff @(posedge clk_i) begin
    if (req_i) begin
      addr_q <= addr_i[AddrBits-1+3:3];
    end
  end

  // this prevents spurious Xes from propagating into
  // the speculative fetch stage of the core
  assign rdata_o = (addr_q < RomSize) ? mem[addr_q] : '0;
endmodule
